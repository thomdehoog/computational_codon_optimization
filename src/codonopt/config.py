"""Design constraints, penalty weights and optimizer settings."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math

from .motifs import Motif, canonicalize_motifs
from .sequence import CODON_TO_AA, validate_dna
from .splice import SpliceScorer
from .tables import CodonUsage, HUMAN

__all__ = ["Weights", "DesignConfig", "OptimizerSettings"]


@dataclass(frozen=True)
class Weights:
    """Relative importance of each penalty term.

    All terms are penalties: lower is better, and zero switches a term off.
    Mandatory constraints are configured separately and take precedence over
    this score. Finite weights alone cannot guarantee a requirement.

    Attributes:
        gc_global: Multiplies absolute CDS GC base-count deviation divided by
            300 (the previous fraction-based scale for a 100-codon CDS).
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
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"weight {name} must be finite and non-negative: {value}")


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
    # Names of mandatory criteria. Other criteria remain soft preferences.
    required: tuple[str, ...] = ("motif",)
    splice_strand: str = "forward"
    left_flank: str = ""
    right_flank: str = ""
    protected_codons: tuple[int, ...] = ()
    excluded_codons: tuple[str, ...] = ()

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
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        for name in ("hairpin_min_loop", "hairpin_max_loop"):
            if isinstance(getattr(self, name), bool) or not isinstance(getattr(self, name), int):
                raise ValueError(f"{name} must be an integer")
        if self.hairpin_min_loop < 0 or self.hairpin_max_loop < self.hairpin_min_loop:
            raise ValueError("require 0 <= hairpin_min_loop <= hairpin_max_loop")
        for name in ("donor_cutoff", "acceptor_cutoff"):
            value = getattr(self, name)
            if value is not None and not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        known = {"motif", "gc_window", "repeat", "hairpin", "homopolymer", "splice"}
        if set(self.required) - known:
            raise ValueError(f"unknown required criteria: {sorted(set(self.required) - known)}")
        object.__setattr__(self, "required", tuple(sorted(set(self.required))))
        if "splice" in self.required and self.splice_scorer is None:
            raise ValueError("required splice screening needs an explicit splice scorer")
        if not isinstance(self.splice_strand, str) or self.splice_strand not in {"forward", "reverse", "both"}:
            raise ValueError("splice_strand must be forward, reverse or both")
        for name in ("left_flank", "right_flank"):
            value = getattr(self, name)
            if not isinstance(value, str):
                raise ValueError(f"{name} must be a DNA string")
            object.__setattr__(self, name, validate_dna(value) if value else "")
        if any(isinstance(i, bool) or not isinstance(i, int) or i < 0 for i in self.protected_codons):
            raise ValueError("protected_codons must contain non-negative integer indices")
        object.__setattr__(self, "protected_codons", tuple(sorted(set(self.protected_codons))))
        excluded = tuple(sorted(set(c.strip().upper().replace("U", "T") for c in self.excluded_codons)))
        if any(c not in CODON_TO_AA for c in excluded):
            raise ValueError("excluded_codons must contain standard genetic-code codons")
        object.__setattr__(self, "excluded_codons", excluded)
        object.__setattr__(self, "motifs", canonicalize_motifs(self.motifs))

    @property
    def splice_strands(self) -> tuple[str, ...]:
        return ("forward", "reverse") if self.splice_strand == "both" else (self.splice_strand,)

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
        should_stop: Called regularly during the search; once it returns True
            the search stops cleanly and returns the best sequence found so
            far, which is how a caller cancels a long run.
        repair_candidates: Maximum complete coordinated candidates per restart.
            Zero disables repair; unresolved requirements are still reported.
        repair_max_codons: Maximum codons per coordinated repair, 2 through 4.
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
    should_stop: object | None = None
    repair_candidates: int = 2000
    repair_max_codons: int = 3

    def __post_init__(self) -> None:
        for name in ("restarts", "max_sweeps", "anneal_steps", "repair_candidates", "repair_max_codons"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f"{name} must be an integer")
        if self.restarts < 1:
            raise ValueError("restarts must be at least 1")
        if self.max_sweeps < 0 or self.anneal_steps < 0:
            raise ValueError("max_sweeps and anneal_steps must not be negative")
        if self.seed_mode not in {"weighted", "best", "random"}:
            raise ValueError(
                f"unknown seed_mode {self.seed_mode!r}; "
                "expected weighted, best or random"
            )
        if self.repair_candidates < 0 or not 2 <= self.repair_max_codons <= 4:
            raise ValueError("repair_candidates must be non-negative and repair_max_codons between 2 and 4")
        if any(not math.isfinite(t) or t <= 0 for t in (self.start_temperature, self.end_temperature)):
            raise ValueError("temperatures must be finite and positive")
        if self.max_seconds is not None and (not math.isfinite(self.max_seconds) or self.max_seconds <= 0):
            raise ValueError("max_seconds must be finite and positive when set")
        if self.seed is not None and (isinstance(self.seed, bool) or not isinstance(self.seed, int)):
            raise ValueError("seed must be an integer")
