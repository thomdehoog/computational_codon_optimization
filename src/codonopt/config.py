"""Design constraints, penalty weights and optimizer settings."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from .motifs import Motif
from .splice import SpliceScorer
from .tables import CodonUsage, HUMAN

__all__ = ["Weights", "DesignConfig", "OptimizerSettings"]


@dataclass(frozen=True)
class Weights:
    """Relative importance of each penalty term.

    All terms are penalties: lower is better, and zero switches a term off.
    The defaults were chosen so that a single forbidden motif outweighs any
    amount of ordinary GC or codon-usage drift, while repeats and hairpins
    trade off against each other roughly one for one.

    Attributes:
        gc_global: Multiplies ``|GC fraction - target|``.
        gc_window: Multiplies the summed shortfall/excess of each scored
            window relative to the acceptable GC range, as a fraction.
        repeat: Multiplies the number of redundant k-mer occurrences.
        hairpin: Multiplies the number of reverse-complement pairs within the
            loop range.
        homopolymer: Multiplies the number of bases by which single-base runs
            exceed the limit.
        motif: Multiplies the number of forbidden motif occurrences.
        splice: Multiplies the summed excess of splice scores over their
            cutoffs.
        usage: Multiplies ``-ln`` of each codon's relative adaptiveness,
            summed over codons.
    """

    gc_global: float = 50.0
    gc_window: float = 5.0
    repeat: float = 1.0
    hairpin: float = 1.0
    homopolymer: float = 2.0
    motif: float = 50.0
    splice: float = 0.5
    usage: float = 0.1

    def __post_init__(self) -> None:
        for name, value in vars(self).items():
            if value < 0:
                raise ValueError(f"weight {name} must not be negative: {value}")


@dataclass(frozen=True)
class DesignConfig:
    """What a good sequence looks like.

    Attributes:
        table: Codon usage table of the expression host.
        gc_target: Overall GC fraction to aim for.
        gc_range: Acceptable GC fraction per window.
        gc_window: Window size in bases for local GC.
        gc_step: Spacing between scored GC windows.
        repeat_k: K-mer size treated as a direct repeat.
        hairpin_stem: Stem length treated as a hairpin arm.
        hairpin_min_loop: Smallest loop counted as a hairpin.
        hairpin_max_loop: Largest loop counted as a hairpin.
        max_homopolymer: Longest acceptable single-base run.
        motifs: Motifs to keep out of the sequence, both strands.
        splice_scorer: Splice backend, or None to skip splice scoring.
        donor_cutoff: Highest acceptable donor score; windows scoring above it
            are penalised. None uses the backend default.
        acceptor_cutoff: Highest acceptable acceptor score. None uses the
            backend default.
        weights: Penalty weights.
    """

    table: CodonUsage = HUMAN
    gc_target: float = 0.50
    gc_range: tuple[float, float] = (0.30, 0.70)
    gc_window: int = 30
    gc_step: int = 3
    repeat_k: int = 12
    hairpin_stem: int = 10
    hairpin_min_loop: int = 3
    hairpin_max_loop: int = 60
    max_homopolymer: int = 6
    motifs: tuple[Motif, ...] = ()
    splice_scorer: SpliceScorer | None = None
    donor_cutoff: float | None = None
    acceptor_cutoff: float | None = None
    weights: Weights = field(default_factory=Weights)

    def __post_init__(self) -> None:
        if not 0.0 <= self.gc_target <= 1.0:
            raise ValueError(
                f"gc_target must be a fraction between 0 and 1: {self.gc_target}"
            )
        low, high = self.gc_range
        if not (0.0 <= low <= high <= 1.0):
            raise ValueError(
                f"gc_range {self.gc_range!r} must be fractions between 0 and 1 "
                "with low <= high (use 0.3 rather than 30)"
            )
        for name in ("gc_window", "gc_step", "repeat_k", "hairpin_stem", "max_homopolymer"):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be at least 1")
        if self.hairpin_min_loop < 0 or self.hairpin_max_loop < self.hairpin_min_loop:
            raise ValueError("require 0 <= hairpin_min_loop <= hairpin_max_loop")

    @property
    def effective_donor_cutoff(self) -> float:
        """Donor cutoff, resolved against the backend default."""
        if self.splice_scorer is None:
            return 0.0
        if self.donor_cutoff is None:
            return self.splice_scorer.default_donor_cutoff
        return self.donor_cutoff

    @property
    def effective_acceptor_cutoff(self) -> float:
        """Acceptor cutoff, resolved against the backend default."""
        if self.splice_scorer is None:
            return 0.0
        if self.acceptor_cutoff is None:
            return self.splice_scorer.default_acceptor_cutoff
        return self.acceptor_cutoff

    def with_(self, **changes) -> "DesignConfig":
        """Return a copy with ``changes`` applied."""
        return replace(self, **changes)


@dataclass(frozen=True)
class OptimizerSettings:
    """How hard to search.

    Attributes:
        restarts: Number of independent seeds; the best result is kept.
        max_sweeps: Greedy passes over every codon position per restart.
        anneal_steps: Simulated-annealing steps run after the greedy passes,
            which lets the search leave a local minimum. Zero skips annealing.
        start_temperature: Initial annealing temperature in score units.
        end_temperature: Final annealing temperature.
        max_seconds: Wall-clock budget per call; the search stops cleanly when
            it is spent. None means no limit.
        seed: Random seed, for reproducible runs.
        seed_mode: How the starting sequence is built -- ``"weighted"`` samples
            codons by host frequency, ``"best"`` takes the most frequent codon
            everywhere, ``"random"`` samples uniformly.
        progress: Called with a short status line after each phase.
    """

    restarts: int = 1
    max_sweeps: int = 6
    anneal_steps: int = 0
    start_temperature: float = 2.0
    end_temperature: float = 0.05
    max_seconds: float | None = 120.0
    seed: int | None = None
    seed_mode: str = "weighted"
    progress: object | None = None

    def __post_init__(self) -> None:
        if self.restarts < 1:
            raise ValueError("restarts must be at least 1")
        if self.max_sweeps < 0 or self.anneal_steps < 0:
            raise ValueError("max_sweeps and anneal_steps must not be negative")
        if self.seed_mode not in {"weighted", "best", "random"}:
            raise ValueError(
                f"unknown seed_mode {self.seed_mode!r}; "
                "expected weighted, best or random"
            )
        if self.start_temperature <= 0 or self.end_temperature <= 0:
            raise ValueError("temperatures must be positive")
        if self.max_seconds is not None and self.max_seconds <= 0:
            raise ValueError("max_seconds must be positive when set")
