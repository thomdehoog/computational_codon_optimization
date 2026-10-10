"""The optimizer.

The original pipeline ran four independent passes in sequence: pick random
codons, fix GC content, remove repeats, then remove splice sites. Each pass
could undo the previous one -- a replacement chosen to break up a repeat was
free to push a window's GC content back out of range, and nothing measured the
overall result. This module replaces the four passes with one weighted
objective (:mod:`codonopt.state`) and a local search over synonymous codons.

Two phases run per restart:

* greedy sweeps, which visit every codon in random order and keep the best
  synonymous replacement, repeating until a whole sweep changes nothing;
* optional simulated annealing, which accepts occasional uphill moves so the
  search can leave a local minimum.

The protein is invariant by construction: the only move available is a swap
between codons of the same amino acid.
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass, field
from typing import Callable

from .config import DesignConfig, OptimizerSettings
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
        initial_sequence: The starting sequence of the restart that produced
            this result, before local search ran on it. With several restarts
            this is the winning restart's seed, not the best seed.
        initial_score: Its objective value.
        initial_breakdown: Its per-term objective values.
        sweeps: Greedy sweeps run, summed over restarts.
        moves_evaluated: Candidate codon swaps evaluated.
        moves_accepted: Swaps that were kept.
        restarts_run: Restarts actually run before the time budget ran out.
        elapsed: Wall-clock seconds spent.
        hit_time_limit: True if the time budget stopped the search early.
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
    notes: list[str] = field(default_factory=list)

    @property
    def improvement(self) -> float:
        """How much the objective fell from the starting sequence."""
        return self.initial_score - self.score


def seed_sequence(
    protein: str,
    table: CodonUsage,
    rng: random.Random,
    mode: str = "weighted",
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

    rng = random.Random(settings.seed)
    started = time.perf_counter()
    deadline = (
        started + settings.max_seconds if settings.max_seconds is not None else None
    )

    best: SequenceState | None = None
    best_initial: tuple[str, float] | None = None
    sweeps = evaluated = accepted = restarts_run = 0
    hit_limit = False

    for restart in range(settings.restarts):
        if deadline is not None and time.perf_counter() >= deadline and best is not None:
            hit_limit = True
            break
        seed = (
            start
            if restart == 0 and start is not None
            else seed_sequence(protein, config.table, rng, settings.seed_mode)
        )
        state = SequenceState(seed, config)
        initial = (state.sequence, state.score)
        restarts_run += 1
        report(
            f"restart {restart + 1}/{settings.restarts}: "
            f"start score {initial[1]:.2f}"
        )

        local_sweeps, local_evaluated, local_accepted, stopped = _greedy(
            state, rng, settings, deadline
        )
        sweeps += local_sweeps
        evaluated += local_evaluated
        accepted += local_accepted
        hit_limit = hit_limit or stopped

        if settings.anneal_steps and not stopped:
            annealed, anneal_evaluated, anneal_accepted, stopped = _anneal(
                state, rng, settings, deadline
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
        if best is None or state.score < best.score:
            best = state
            best_initial = initial

    assert best is not None and best_initial is not None  # restarts >= 1
    if translate(best.sequence) != protein:  # pragma: no cover - safety net
        raise AssertionError("optimizer changed the encoded protein")

    initial_state = SequenceState(best_initial[0], config)
    return OptimizationResult(
        sequence=best.sequence,
        protein=protein,
        score=best.score,
        breakdown=best.breakdown(),
        initial_sequence=best_initial[0],
        initial_score=best_initial[1],
        initial_breakdown=initial_state.breakdown(),
        sweeps=sweeps,
        moves_evaluated=evaluated,
        moves_accepted=accepted,
        restarts_run=restarts_run,
        elapsed=time.perf_counter() - started,
        hit_time_limit=hit_limit,
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
    deadline: float | None,
) -> tuple[int, int, int, bool]:
    """Run greedy sweeps until nothing improves. Returns counters."""
    order = list(range(state.codon_count))
    sweeps = evaluated = accepted = 0
    for _ in range(settings.max_sweeps):
        rng.shuffle(order)
        sweeps += 1
        improved = 0
        for counter, position in enumerate(order):
            if deadline is not None and not counter % 64:
                if time.perf_counter() >= deadline:
                    return sweeps, evaluated, accepted, True
            options = state.alternatives(position)
            if len(options) == 1:
                continue
            current = state.codon(position)
            best_codon, best_score = current, state.score
            for option in options:
                if option == current:
                    continue
                state.set_codon(position, option)
                evaluated += 1
                score = state.score
                if score < best_score - _EPSILON:
                    best_codon, best_score = option, score
            state.set_codon(position, best_codon)
            if best_codon != current:
                improved += 1
        accepted += improved
        if not improved:
            break
    return sweeps, evaluated, accepted, False


def _anneal(
    state: SequenceState,
    rng: random.Random,
    settings: OptimizerSettings,
    deadline: float | None,
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

    best_sequence, best_score = state.sequence, state.score
    cooling = (settings.end_temperature / settings.start_temperature) ** (
        1.0 / max(1, steps - 1)
    )
    temperature = settings.start_temperature
    evaluated = accepted = 0
    stopped = False

    for step in range(steps):
        if deadline is not None and not step % 256 and time.perf_counter() >= deadline:
            stopped = True
            break
        position = rng.choice(positions)
        current = state.codon(position)
        options = state.alternatives(position)
        option = rng.choice(options)
        if option == current:
            continue
        previous_score = state.score
        state.set_codon(position, option)
        evaluated += 1
        delta = state.score - previous_score
        if delta <= 0 or rng.random() < math.exp(-delta / temperature):
            accepted += 1
            if state.score < best_score - _EPSILON:
                best_sequence, best_score = state.sequence, state.score
        else:
            state.set_codon(position, current)
        temperature *= cooling

    if state.sequence != best_sequence:
        state = SequenceState(best_sequence, state.config)
    return state, evaluated, accepted, stopped
