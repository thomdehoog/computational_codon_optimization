"""The optimizer.

The original pipeline ran four independent passes in sequence: pick random
codons, fix GC content, remove repeats, then remove splice sites. Each pass
could undo the previous one -- a replacement chosen to break up a repeat was
free to push a window's GC content back out of range, and nothing measured the
overall result. This module replaces the four passes with one weighted
objective (:mod:`codonopt.state`) and a local search over synonymous codons.

Three phases run per restart:

* greedy sweeps, which visit every codon in random order and keep the best
  synonymous replacement, repeating until a whole sweep changes nothing;
* bounded coordinated repair around unresolved mandatory requirements;
* optional simulated annealing, which accepts occasional uphill moves so the
  search can leave a local minimum.

The protein is invariant by construction: the only move available is a swap
between codons of the same amino acid.
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass, field, fields
from itertools import product
from typing import Callable

from .config import DesignConfig, OptimizerSettings
from .constraints import constraint_violations
from .sequence import (
    AA_TO_CODONS,
    translate,
    validate_dna,
    validate_protein,
)
from .state import ScoreBreakdown, SequenceState
from .tables import CodonUsage

__all__ = ["OptimizationResult", "seed_sequence", "optimize", "optimize_dna"]

_EPSILON = 1e-12


@dataclass
class OptimizationResult:
    """What a run produced.

    Attributes:
        sequence: The optimized coding sequence.
        protein: The protein it encodes, which equals the input protein.
        score: Final objective value.
        breakdown: Final per-term objective values.
        initial_sequence: Original DNA input, or the first seed for a protein
            design. This baseline is preserved even if a later restart wins.
        initial_score: Its objective value.
        initial_breakdown: Its per-term objective values.
        sweeps: Greedy sweeps run, summed over restarts.
        moves_evaluated: Candidate codon swaps evaluated.
        moves_accepted: Swaps that were kept.
        restarts_run: Restarts actually run before the time budget ran out.
        elapsed: Wall-clock seconds spent.
        hit_time_limit: True if the time budget stopped the search early.
        cancelled: True if ``should_stop`` stopped the search early.
    """

    sequence: str
    protein: str
    score: float
    breakdown: ScoreBreakdown
    initial_sequence: str
    initial_score: float
    initial_breakdown: ScoreBreakdown
    sweeps: int = 0
    moves_evaluated: int = 0
    moves_accepted: int = 0
    restarts_run: int = 0
    elapsed: float = 0.0
    hit_time_limit: bool = False
    cancelled: bool = False
    notes: list[str] = field(default_factory=list)
    valid: bool = True
    unresolved_constraints: list[dict] = field(default_factory=list)
    stop_reason: str = "completed"
    effective_seed: int | None = None
    search_settings: dict = field(default_factory=dict)
    winning_seed_sequence: str = ""
    repair_evaluated: int = 0

    @property
    def improvement(self) -> float:
        """How much the objective fell from the starting sequence."""
        return self.initial_score - self.score


def seed_sequence(
    protein: str,
    table: CodonUsage,
    rng: random.Random,
    mode: str = "weighted",
    *,
    excluded_codons: tuple[str, ...] = (),
) -> str:
    """Build a starting coding sequence for ``protein``.

    Args:
        protein: Protein sequence, optionally ending in ``*``.
        table: Codon usage table of the host.
        rng: Random source.
        mode: ``"weighted"`` samples codons in proportion to host frequency,
            ``"best"`` always takes the most frequent codon, ``"random"``
            samples uniformly as the original scripts did.

    Raises:
        ValueError: On an unknown ``mode``.
    """
    if excluded_codons:
        sequence = []
        for aa in protein:
            options = [c for c in AA_TO_CODONS[aa] if c not in excluded_codons]
            if not options:
                raise ValueError(f"all codons for {aa} are excluded")
            if mode == "weighted":
                frequencies = [table.frequency(c) for c in options]
                codon = rng.choices(options, weights=frequencies if sum(frequencies) else None, k=1)[0]
            elif mode == "best":
                codon = max(options, key=table.frequency)
            elif mode == "random":
                codon = rng.choice(options)
            else:
                raise ValueError(f"unknown seed mode {mode!r}")
            sequence.append(codon)
        return "".join(sequence)
    if mode == "weighted":
        return "".join(table.sample(aa, rng) for aa in protein)
    if mode == "best":
        return "".join(table.most_frequent(aa) for aa in protein)
    if mode == "random":
        return "".join(rng.choice(AA_TO_CODONS[aa]) for aa in protein)
    raise ValueError(f"unknown seed mode {mode!r}")


def optimize(
    protein: str,
    config: DesignConfig | None = None,
    settings: OptimizerSettings | None = None,
    *,
    start: str | None = None,
) -> OptimizationResult:
    """Design a coding sequence for ``protein``.

    Args:
        protein: Protein sequence, optionally ending in ``*``.
        config: Design constraints and weights; defaults apply if omitted.
        settings: Search effort; defaults apply if omitted.
        start: Optional in-frame coding sequence to start the first restart
            from, for example an existing construct to improve. It must encode
            ``protein``.

    Raises:
        codonopt.sequence.SequenceError: If ``protein`` or ``start`` is invalid.
        ValueError: If ``start`` does not encode ``protein``.
    """
    protein = validate_protein(protein)
    config = config or DesignConfig()
    settings = settings or OptimizerSettings()
    report: Callable[[str], None] = settings.progress or (lambda message: None)  # type: ignore[assignment]

    if start is not None:
        start = validate_dna(start, in_frame=True)
        if translate(start) != protein:
            raise ValueError("start sequence does not encode the given protein")

    if config.protected_codons and start is None:
        raise ValueError("protected_codons require a starting DNA sequence")
    if any(i >= len(protein) for i in config.protected_codons):
        raise ValueError("protected codon index is outside the coding sequence")
    for i, aa in enumerate(protein):
        if i not in config.protected_codons and not set(AA_TO_CODONS[aa]) - set(config.excluded_codons):
            raise ValueError(f"all codons for {aa} at position {i} are excluded")

    effective_seed = settings.seed if settings.seed is not None else random.SystemRandom().getrandbits(31)
    rng = random.Random(effective_seed)
    started = time.perf_counter()
    deadline = (
        started + settings.max_seconds if settings.max_seconds is not None else None
    )
    should_stop: Callable[[], bool] | None = settings.should_stop  # type: ignore[assignment]
    stop_reason = "completed"

    def expired() -> bool:
        """True once the time budget is spent or the caller asked to stop."""
        nonlocal stop_reason
        if stop_reason != "completed":
            return True
        if should_stop is not None and should_stop():
            stop_reason = "cancelled"
            return True
        if deadline is not None and time.perf_counter() >= deadline:
            stop_reason = "time_limit"
            return True
        return False

    best: SequenceState | None = None
    winning_seed = ""
    baseline: SequenceState | None = None
    sweeps = evaluated = accepted = restarts_run = 0
    hit_limit = False
    repair_evaluated = 0

    for restart in range(settings.restarts):
        if best is not None and expired():
            hit_limit = True
            break
        if restart == 0 and start is not None:
            seed = start
        elif config.protected_codons:
            seed = "".join(
                start[3*i:3*i+3] if i in config.protected_codons else
                seed_sequence(aa, config.table, rng, settings.seed_mode,
                              excluded_codons=config.excluded_codons)
                for i, aa in enumerate(protein)
            )
        else:
            seed = seed_sequence(protein, config.table, rng, settings.seed_mode,
                                 excluded_codons=config.excluded_codons)
        if config.protected_codons:
            codons = [seed[i:i+3] for i in range(0, len(seed), 3)]
            for i in config.protected_codons:
                codons[i] = start[3*i:3*i+3]
            seed = "".join(codons)
        state = SequenceState(seed, config)
        if baseline is None:
            baseline = SequenceState(seed, config)
        initial = (state.sequence, state.score)
        restarts_run += 1
        report(
            f"restart {restart + 1}/{settings.restarts}: "
            f"start score {initial[1]:.2f}"
        )

        local_sweeps, local_evaluated, local_accepted, stopped = _greedy(
            state, rng, settings, expired
        )
        sweeps += local_sweeps
        evaluated += local_evaluated
        accepted += local_accepted
        hit_limit = hit_limit or stopped

        if state.violations and not stopped and settings.repair_candidates:
            count, kept, stopped = _repair(state, settings, expired)
            repair_evaluated += count
            evaluated += count
            accepted += kept
            hit_limit = hit_limit or stopped
            report(f"repair: {count} candidates, {state.violations} mandatory violations remain")
            if kept and not stopped:
                extra_sweeps, extra_evaluated, extra_accepted, stopped = _greedy(state, rng, settings, expired)
                sweeps += extra_sweeps
                evaluated += extra_evaluated
                accepted += extra_accepted
                hit_limit = hit_limit or stopped

        if settings.anneal_steps and not stopped:
            annealed, anneal_evaluated, anneal_accepted, stopped = _anneal(
                state, rng, settings, expired
            )
            evaluated += anneal_evaluated
            accepted += anneal_accepted
            hit_limit = hit_limit or stopped
            if annealed is not None:
                state = annealed

        report(
            f"restart {restart + 1}/{settings.restarts}: "
            f"final score {state.score:.2f}"
        )
        if best is None or _better(state.rank, best.rank):
            best = state
            winning_seed = seed

    assert best is not None and baseline is not None  # restarts >= 1
    if translate(best.sequence) != protein:  # pragma: no cover - safety net
        raise AssertionError("optimizer changed the encoded protein")

    # Validate independently rather than treating the incremental counter as
    # proof of feasibility. Never label heuristic failure as impossibility.
    unresolved = constraint_violations(best.sequence, config)
    if bool(unresolved) != bool(best.violations):
        raise AssertionError("incremental constraint count disagrees with independent validation")
    cancelled = stop_reason == "cancelled"
    search_settings = {f.name: getattr(settings, f.name) for f in fields(settings)
                       if f.name not in {"progress", "should_stop"}}
    search_settings["seed"] = effective_seed
    return OptimizationResult(
        sequence=best.sequence,
        protein=protein,
        score=best.score,
        breakdown=best.breakdown(),
        initial_sequence=baseline.sequence,
        initial_score=baseline.score,
        initial_breakdown=baseline.breakdown(),
        sweeps=sweeps,
        moves_evaluated=evaluated,
        moves_accepted=accepted,
        restarts_run=restarts_run,
        elapsed=time.perf_counter() - started,
        hit_time_limit=hit_limit and not cancelled,
        cancelled=cancelled,
        valid=not unresolved,
        unresolved_constraints=unresolved,
        stop_reason=stop_reason,
        effective_seed=effective_seed,
        search_settings=search_settings,
        winning_seed_sequence=winning_seed,
        repair_evaluated=repair_evaluated,
        notes=(["Mandatory constraints remain unresolved; this is not proof of impossibility."] if unresolved else []),
    )


def optimize_dna(
    sequence: str,
    config: DesignConfig | None = None,
    settings: OptimizerSettings | None = None,
) -> OptimizationResult:
    """Improve an existing coding sequence without changing its protein."""
    sequence = validate_dna(sequence, in_frame=True)
    return optimize(
        translate(sequence), config, settings, start=sequence
    )


def _greedy(
    state: SequenceState,
    rng: random.Random,
    settings: OptimizerSettings,
    expired: Callable[[], bool],
) -> tuple[int, int, int, bool]:
    """Run greedy sweeps until nothing improves. Returns counters."""
    if expired():
        return 0, 0, 0, True
    order = list(range(state.codon_count))
    sweeps = evaluated = accepted = 0
    for _ in range(settings.max_sweeps):
        rng.shuffle(order)
        sweeps += 1
        improved = 0
        for counter, position in enumerate(order):
            if not counter % 64 and expired():
                return sweeps, evaluated, accepted + improved, True
            options = state.alternatives(position)
            if options == (state.codon(position),):
                continue
            current = state.codon(position)
            best_codon, best_rank = current, state.rank
            for option in options:
                if option == current:
                    continue
                state.set_codon(position, option)
                evaluated += 1
                rank = state.rank
                if _better(rank, best_rank):
                    best_codon, best_rank = option, rank
            state.set_codon(position, best_codon)
            if best_codon != current:
                improved += 1
        accepted += improved
        if not improved:
            break
    return sweeps, evaluated, accepted, False


def _better(candidate: tuple[int, float], incumbent: tuple[int, float]) -> bool:
    return candidate[0] < incumbent[0] or (candidate[0] == incumbent[0] and candidate[1] < incumbent[1] - _EPSILON)


def _repair(state: SequenceState, settings: OptimizerSettings,
            expired: Callable[[], bool]) -> tuple[int, int, bool]:
    """Bounded coordinated repair around unresolved regions.

    Intermediate edits may violate constraints; only the complete candidate
    is compared with the incumbent. The search budget is per restart and is
    shared across passes and regions.
    """
    evaluated = accepted = 0
    while state.violations and evaluated < settings.repair_candidates:
        changed = False
        failures = constraint_violations(state.sequence, state.config)
        for failure in failures:
            for region in failure["regions"]:
                offset = len(state.config.left_flank)
                lo = max(0, (region["start"] - offset) // 3)
                hi = min(state.codon_count, (region["end"] - offset + 2) // 3)
                if lo >= hi:  # wholly fixed context cannot be edited
                    continue
                mutable = [p for p in range(max(0, lo - 1), min(state.codon_count, hi + 1))
                           if len(state.alternatives(p)) > 1]
                # Long failures are addressed with overlapping local blocks;
                # do not enumerate combinations across an entire long run.
                for width in range(2, settings.repair_max_codons + 1):
                    for block_start in range(max(1, len(mutable) - width + 1)):
                        positions = mutable[block_start:block_start + width]
                        if len(positions) != width:
                            continue
                        original = tuple(state.codon(p) for p in positions)
                        best_codons, best_rank = original, state.rank
                        options = [state.alternatives(p) for p in positions]
                        for candidate in product(*options):
                            if candidate == original:
                                continue
                            if evaluated >= settings.repair_candidates or expired():
                                for p, c in zip(positions, best_codons):
                                    state.set_codon(p, c)
                                return evaluated, accepted + int(best_codons != original), expired()
                            for p, c in zip(positions, candidate):
                                state.set_codon(p, c)
                            evaluated += 1
                            if _better(state.rank, best_rank):
                                best_codons, best_rank = candidate, state.rank
                        for p, c in zip(positions, best_codons):
                            state.set_codon(p, c)
                        if best_codons != original:
                            accepted += 1
                            changed = True
                        if state.violations == 0:
                            return evaluated, accepted, False
        if not changed:
            break
    return evaluated, accepted, False


def _anneal(
    state: SequenceState,
    rng: random.Random,
    settings: OptimizerSettings,
    expired: Callable[[], bool],
) -> tuple[SequenceState | None, int, int, bool]:
    """Run simulated annealing, returning the best state seen.

    The incumbent is allowed to drift uphill, so the best sequence seen is
    tracked separately and rebuilt at the end.
    """
    steps = settings.anneal_steps
    positions = [
        index
        for index in range(state.codon_count)
        if len(state.alternatives(index)) > 1
    ]
    if not positions or steps <= 0:
        return None, 0, 0, False

    best_sequence, best_rank = state.sequence, state.rank
    evaluated = accepted = 0
    stopped = False

    for step in range(steps):
        # Schedule is tied to iterations, including draws of the current codon.
        fraction = step / max(1, steps - 1)
        temperature = math.exp((1-fraction)*math.log(settings.start_temperature)
                               + fraction*math.log(settings.end_temperature))
        if not step % 256 and expired():
            stopped = True
            break
        position = rng.choice(positions)
        current = state.codon(position)
        options = state.alternatives(position)
        option = rng.choice(options)
        if option == current:
            continue
        previous_rank = state.rank
        state.set_codon(position, option)
        evaluated += 1
        delta = state.score - previous_rank[1]
        if state.violations < previous_rank[0] or (
            state.violations == previous_rank[0] and
            (delta <= 0 or rng.random() < math.exp(-delta / temperature))
        ):
            accepted += 1
            if _better(state.rank, best_rank):
                best_sequence, best_rank = state.sequence, state.rank
        else:
            state.set_codon(position, current)

    if state.sequence != best_sequence:
        state = SequenceState(best_sequence, state.config)
    return state, evaluated, accepted, stopped
