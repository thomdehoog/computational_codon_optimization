"""Incrementally scored sequence state.

:class:`SequenceState` holds a coding sequence together with the value of the
weighted penalty objective. Every penalty term is attributed to a position --
a GC window start, a base that extends a homopolymer run, a motif match, a
splice window, a codon -- so after a codon swap only the terms anchored near
the edit need recomputing. Repeats and hairpins, which are not local, are
carried by :class:`codonopt.index.KmerIndex`.

The result is that a candidate codon swap costs time proportional to the window
sizes involved rather than to the length of the sequence, which is what makes a
whole-sequence search practical. :func:`SequenceState.audit` recomputes the
score from scratch and is used by the tests to confirm the incremental path
agrees with a full recomputation.
"""

from __future__ import annotations

import math
import re
from bisect import bisect_left, bisect_right
from dataclasses import dataclass

from .config import DesignConfig
from .index import KmerIndex
from .metrics import window_starts
from .motifs import Motif
from .sequence import (
    AA_TO_CODONS,
    CODON_TO_AA,
    gc_count,
    reverse_complement,
    split_codons,
    translate,
    validate_dna,
)
from .splice import ACCEPTOR_WINDOW, DONOR_WINDOW

__all__ = ["SequenceState", "ScoreBreakdown"]

_MIN_WEIGHT = 1e-9


@dataclass(frozen=True)
class ScoreBreakdown:
    """Per-term view of the objective, for reports.

    Attributes:
        total: The weighted objective value.
        terms: Term name -> weighted contribution.
        raw: Term name -> the unweighted quantity being penalised.
    """

    total: float
    terms: dict[str, float]
    raw: dict[str, float]


class SequenceState:
    """A coding sequence plus its incrementally maintained penalty score.

    Args:
        sequence: In-frame coding sequence.
        config: Design constraints and weights.

    Raises:
        codonopt.sequence.SequenceError: If ``sequence`` is not in-frame DNA.
    """

    __slots__ = (
        "config",
        "protein",
        "_sequence",
        "_gc",
        "_local_total",
        "_repeat_index",
        "_hairpin_index",
        "_window_starts",
        "_usage_cost",
        "_literals_by_length",
        "_pattern_motifs",
        "_motif_lengths",
        "_splice_terms",
    )

    def __init__(self, sequence: str, config: DesignConfig) -> None:
        self._sequence = validate_dna(sequence, in_frame=True)
        self.config = config
        self.protein = translate(self._sequence)
        weights = config.weights

        self._gc = gc_count(self._sequence)
        self._window_starts = window_starts(
            len(self._sequence), config.gc_window, config.gc_step
        )
        self._usage_cost = {
            codon: -math.log(max(config.table.weight(codon), 1e-6))
            for codon in CODON_TO_AA
        }
        self._literals_by_length, self._pattern_motifs = _split_motifs(config.motifs)
        self._motif_lengths = sorted(self._literals_by_length)
        self._splice_terms = (
            (
                (DONOR_WINDOW, config.splice_scorer.score_donor,
                 config.effective_donor_cutoff),
                (ACCEPTOR_WINDOW, config.splice_scorer.score_acceptor,
                 config.effective_acceptor_cutoff),
            )
            if config.splice_scorer is not None and weights.splice > _MIN_WEIGHT
            else ()
        )

        self._repeat_index = (
            KmerIndex(self._sequence, config.repeat_k, track_hairpins=False)
            if weights.repeat > _MIN_WEIGHT
            else None
        )
        self._hairpin_index = (
            KmerIndex(
                self._sequence,
                config.hairpin_stem,
                min_loop=config.hairpin_min_loop,
                max_loop=config.hairpin_max_loop,
                track_hairpins=True,
            )
            if weights.hairpin > _MIN_WEIGHT
            else None
        )
        self._local_total = self._local(0, len(self._sequence))

    # -- read-only views ---------------------------------------------------

    @property
    def sequence(self) -> str:
        """The current coding sequence."""
        return self._sequence

    @property
    def codon_count(self) -> int:
        """Number of codons in the sequence."""
        return len(self._sequence) // 3

    @property
    def score(self) -> float:
        """The weighted objective; lower is better."""
        return self._local_total + self._global_terms()

    def codon(self, codon_index: int) -> str:
        """Return the codon at ``codon_index``."""
        start = codon_index * 3
        return self._sequence[start : start + 3]

    def alternatives(self, codon_index: int) -> tuple[str, ...]:
        """Return every synonymous codon for the residue at ``codon_index``."""
        return AA_TO_CODONS[CODON_TO_AA[self.codon(codon_index)]]

    # -- mutation ----------------------------------------------------------

    def set_codon(self, codon_index: int, codon: str) -> str:
        """Replace one codon in place and update the score.

        The replacement must be synonymous, so the encoded protein cannot
        change. Returns the codon that was there before, which is all
        :meth:`set_codon` needs to undo an edit.

        Raises:
            ValueError: If ``codon`` does not encode the same residue.
        """
        lo = codon_index * 3
        hi = lo + 3
        previous = self._sequence[lo:hi]
        if codon == previous:
            return previous
        if CODON_TO_AA[codon] != CODON_TO_AA[previous]:
            raise ValueError(
                f"{codon} encodes {CODON_TO_AA[codon]} but position "
                f"{codon_index} holds {CODON_TO_AA[previous]}"
            )

        before = self._local(lo, hi)
        self._sequence = self._sequence[:lo] + codon + self._sequence[hi:]
        self._gc += gc_count(codon) - gc_count(previous)
        if self._repeat_index is not None:
            self._repeat_index.replace_range(self._sequence, lo, hi)
        if self._hairpin_index is not None:
            self._hairpin_index.replace_range(self._sequence, lo, hi)
        self._local_total += self._local(lo, hi) - before
        return previous

    # -- scoring -----------------------------------------------------------

    def _global_terms(self) -> float:
        weights = self.config.weights
        total = 0.0
        if weights.gc_global > _MIN_WEIGHT:
            fraction = self._gc / len(self._sequence)
            total += weights.gc_global * abs(fraction - self.config.gc_target)
        if self._repeat_index is not None:
            total += weights.repeat * self._repeat_index.redundant
        if self._hairpin_index is not None:
            total += weights.hairpin * self._hairpin_index.hairpin_pairs
        return total

    def _local(self, lo: int, hi: int) -> float:
        """Weighted penalty of every term anchored near bases ``[lo, hi)``."""
        sequence = self._sequence
        length = len(sequence)
        config = self.config
        weights = config.weights
        total = 0.0

        if weights.gc_window > _MIN_WEIGHT and self._window_starts:
            window = config.gc_window
            low_bound, high_bound = config.gc_range
            starts = self._window_starts
            first = bisect_left(starts, lo - window + 1)
            last = bisect_right(starts, hi - 1)
            for position in starts[first:last]:
                end = min(position + window, length)
                fraction = (
                    sequence.count("G", position, end)
                    + sequence.count("C", position, end)
                ) / (end - position)
                if fraction < low_bound:
                    total += weights.gc_window * (low_bound - fraction)
                elif fraction > high_bound:
                    total += weights.gc_window * (fraction - high_bound)

        if weights.homopolymer > _MIN_WEIGHT:
            run = config.max_homopolymer
            for i in range(max(run, lo), min(length, hi + run)):
                segment = sequence[i - run : i + 1]
                if segment == segment[0] * (run + 1):
                    total += weights.homopolymer

        if weights.motif > _MIN_WEIGHT:
            for size in self._motif_lengths:
                literals = self._literals_by_length[size]
                first = max(0, lo - size + 1)
                last = min(length - size + 1, hi)
                for position in range(first, last):
                    if sequence[position : position + size] in literals:
                        total += weights.motif
            for motif in self._pattern_motifs:
                first = max(0, lo - motif.length + 1)
                last = min(length - motif.length + 1, hi)
                for position in range(first, last):
                    if motif.matches_at(sequence, position):
                        total += weights.motif

        for window, score_window, cutoff in self._splice_terms:
            first = max(0, lo - window + 1)
            last = min(length - window + 1, hi)
            weight = weights.splice
            for position in range(first, last):
                forward = sequence[position : position + window]
                score = score_window(forward)
                if score > cutoff:
                    total += weight * (score - cutoff)
                score = score_window(reverse_complement(forward))
                if score > cutoff:
                    total += weight * (score - cutoff)

        if weights.usage > _MIN_WEIGHT:
            cost = self._usage_cost
            for index in range(lo // 3, (hi + 2) // 3):
                total += weights.usage * cost[sequence[index * 3 : index * 3 + 3]]

        return total

    # -- reporting ---------------------------------------------------------

    def audit(self) -> float:
        """Return the score recomputed from scratch.

        Used by the tests to confirm the incremental bookkeeping has not
        drifted from a full recomputation.
        """
        return SequenceState(self._sequence, self.config).score

    def breakdown(self) -> ScoreBreakdown:
        """Return the objective split into its individual terms."""
        config = self.config
        weights = config.weights
        sequence = self._sequence
        raw: dict[str, float] = {}
        terms: dict[str, float] = {}

        fraction = self._gc / len(sequence)
        raw["gc_global"] = fraction
        terms["gc_global"] = weights.gc_global * abs(fraction - config.gc_target)

        window_excess = 0.0
        low_bound, high_bound = config.gc_range
        for position in self._window_starts:
            end = min(position + config.gc_window, len(sequence))
            local = (
                sequence.count("G", position, end) + sequence.count("C", position, end)
            ) / (end - position)
            if local < low_bound:
                window_excess += low_bound - local
            elif local > high_bound:
                window_excess += local - high_bound
        raw["gc_window"] = window_excess
        terms["gc_window"] = weights.gc_window * window_excess

        repeat_index = self._repeat_index or KmerIndex(
            sequence, config.repeat_k, track_hairpins=False
        )
        raw["repeat"] = float(repeat_index.redundant)
        terms["repeat"] = weights.repeat * repeat_index.redundant

        hairpin_index = self._hairpin_index or KmerIndex(
            sequence,
            config.hairpin_stem,
            min_loop=config.hairpin_min_loop,
            max_loop=config.hairpin_max_loop,
            track_hairpins=True,
        )
        raw["hairpin"] = float(hairpin_index.hairpin_pairs)
        terms["hairpin"] = weights.hairpin * hairpin_index.hairpin_pairs

        run = config.max_homopolymer
        excess = sum(
            1
            for i in range(run, len(sequence))
            if sequence[i - run : i + 1] == sequence[i] * (run + 1)
        )
        raw["homopolymer"] = float(excess)
        terms["homopolymer"] = weights.homopolymer * excess

        hits = 0
        for motif in config.motifs:
            hits += len(motif.find_all(sequence))
        raw["motif"] = float(hits)
        terms["motif"] = weights.motif * hits

        splice_excess = 0.0
        for window, score_window, cutoff in self._splice_terms:
            for position in range(len(sequence) - window + 1):
                forward = sequence[position : position + window]
                for candidate in (forward, reverse_complement(forward)):
                    score = score_window(candidate)
                    if score > cutoff:
                        splice_excess += score - cutoff
        raw["splice"] = splice_excess
        terms["splice"] = weights.splice * splice_excess

        usage = sum(self._usage_cost[codon] for codon in split_codons(sequence))
        raw["usage"] = usage
        terms["usage"] = weights.usage * usage

        return ScoreBreakdown(
            total=sum(terms.values()), terms=terms, raw=raw
        )


def _split_motifs(
    motifs: tuple[Motif, ...] | list[Motif],
) -> tuple[dict[int, frozenset[str]], list[Motif]]:
    """Split motifs into literal strings by length and IUPAC patterns.

    Literal motifs, which are the great majority, become set-membership tests
    on fixed-length slices -- much cheaper in the inner loop than running a
    regular expression per motif per position. Both strands are included.
    """
    literals: dict[int, set[str]] = {}
    patterns: list[Motif] = []
    for motif in motifs:
        if re.fullmatch(r"[ACGT]+", motif.pattern):
            for variant in {motif.pattern, reverse_complement(motif.pattern)}:
                literals.setdefault(len(variant), set()).add(variant)
        else:
            patterns.append(motif)
    return {size: frozenset(values) for size, values in literals.items()}, patterns
