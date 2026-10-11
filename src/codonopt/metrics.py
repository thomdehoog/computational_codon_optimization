"""Sequence quality measurements used for reporting and as test oracles.

These functions recompute each property from scratch over a whole sequence.
The optimizer does not call them in its inner loop -- it maintains the same
quantities incrementally (see :mod:`codonopt.state`) -- but they define what
those incremental counters are supposed to mean, and the test suite checks the
two agree.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from bisect import bisect_left, bisect_right

from .motifs import Motif
from .sequence import gc_count, gc_prefix, reverse_complement

__all__ = [
    "GcWindowReport",
    "HomopolymerRun",
    "RepeatReport",
    "Hairpin",
    "gc_profile",
    "gc_window_report",
    "homopolymer_runs",
    "homopolymer_excess",
    "repeat_report",
    "find_hairpins",
    "count_hairpin_pairs",
    "motif_hits",
    "window_starts",
]


def window_starts(length: int, window: int, step: int) -> list[int]:
    """Return the start offsets of every scored window.

    Windows are ``window`` bases long and spaced ``step`` apart. The final
    window is anchored to the end of the sequence so the tail is never left
    unscored, and no window extends past the sequence.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window and step must be positive")
    if length < window:
        return [0] if length else []
    starts = list(range(0, length - window + 1, step))
    if starts[-1] != length - window:
        starts.append(length - window)
    return starts


@dataclass(frozen=True)
class GcWindowReport:
    """Local GC content statistics.

    Attributes:
        overall: GC fraction of the whole sequence.
        window: Window size in bases.
        step: Spacing between scored windows.
        low: Number of scored windows below the target range.
        high: Number of scored windows above the target range.
        windows: Number of scored windows.
        min_fraction: Lowest window GC fraction.
        max_fraction: Highest window GC fraction.
        worst_positions: 0-based starts of the offending windows.
    """

    overall: float
    window: int
    step: int
    low: int
    high: int
    windows: int
    min_fraction: float
    max_fraction: float
    worst_positions: tuple[int, ...] = ()

    @property
    def violations(self) -> int:
        """Number of scored windows outside the target range."""
        return self.low + self.high


def gc_window_report(
    sequence: str,
    *,
    window: int = 30,
    step: int = 3,
    target_range: tuple[float, float] = (0.30, 0.70),
) -> GcWindowReport:
    """Measure overall and windowed GC content.

    Args:
        sequence: Nucleotide sequence.
        window: Window size in bases.
        step: Spacing between scored windows.
        target_range: Acceptable window GC content, as *fractions* between 0
            and 1. Percentages are rejected rather than silently misread,
            which is the bug the original scripts carried.

    Raises:
        ValueError: If ``target_range`` is not a pair of fractions in [0, 1]
            with low <= high.
    """
    low_bound, high_bound = _check_fraction_range(target_range)
    if not sequence:
        return GcWindowReport(0.0, window, step, 0, 0, 0, 0.0, 0.0)

    prefix = gc_prefix(sequence)
    starts = window_starts(len(sequence), window, step)
    low = high = 0
    worst: list[int] = []
    lowest, highest = 1.0, 0.0
    for start in starts:
        end = min(start + window, len(sequence))
        fraction = (prefix[end] - prefix[start]) / (end - start)
        lowest = min(lowest, fraction)
        highest = max(highest, fraction)
        if fraction < low_bound:
            low += 1
            worst.append(start)
        elif fraction > high_bound:
            high += 1
            worst.append(start)
    return GcWindowReport(
        overall=gc_count(sequence) / len(sequence),
        window=window,
        step=step,
        low=low,
        high=high,
        windows=len(starts),
        min_fraction=lowest,
        max_fraction=highest,
        worst_positions=tuple(worst),
    )


def gc_profile(
    sequence: str, *, window: int = 30, step: int = 3
) -> list[tuple[int, float]]:
    """Return ``(start, GC fraction)`` for every scored window.

    The same windows :func:`gc_window_report` scores, so a plot of this profile
    shows exactly what the penalty sees.
    """
    if not sequence:
        return []
    prefix = gc_prefix(sequence)
    profile = []
    for start in window_starts(len(sequence), window, step):
        end = min(start + window, len(sequence))
        profile.append((start, (prefix[end] - prefix[start]) / (end - start)))
    return profile


def _check_fraction_range(target_range: tuple[float, float]) -> tuple[float, float]:
    low, high = target_range
    if not (0.0 <= low <= high <= 1.0):
        raise ValueError(
            f"target GC range {target_range!r} must be fractions between 0 and 1 "
            "with low <= high (use 0.3 rather than 30)"
        )
    return float(low), float(high)


@dataclass(frozen=True)
class HomopolymerRun:
    """A run of one repeated base.

    Attributes:
        base: The repeated base.
        start: 0-based start position.
        length: Run length in bases.
    """

    base: str
    start: int
    length: int


def homopolymer_runs(sequence: str, min_length: int = 7) -> list[HomopolymerRun]:
    """Return every maximal single-base run of at least ``min_length`` bases."""
    if min_length < 1:
        raise ValueError("min_length must be at least 1")
    runs: list[HomopolymerRun] = []
    start = 0
    for i in range(1, len(sequence) + 1):
        if i == len(sequence) or sequence[i] != sequence[start]:
            length = i - start
            if length >= min_length:
                runs.append(HomopolymerRun(sequence[start], start, length))
            start = i
    return runs


def homopolymer_excess(sequence: str, max_run: int = 6) -> int:
    """Return the number of bases by which runs exceed ``max_run``.

    This counts one unit for every position that extends a run past the limit,
    which is the same quantity the optimizer maintains incrementally, and it
    equals ``sum(len - max_run for runs longer than max_run)``.
    """
    if max_run < 1:
        raise ValueError("max_run must be at least 1")
    excess = 0
    for i in range(max_run, len(sequence)):
        window = sequence[i - max_run : i + 1]
        if window == window[0] * (max_run + 1):
            excess += 1
    return excess


@dataclass(frozen=True)
class RepeatReport:
    """Direct-repeat statistics for one k-mer size.

    Attributes:
        k: K-mer size in bases.
        redundant: Number of k-mer occurrences beyond the first for each
            distinct repeated k-mer. Zero means every k-mer is unique.
        distinct_repeated: How many distinct k-mers occur more than once.
        max_occurrences: Highest occurrence count of any k-mer.
        covered_fraction: Fraction of sequence positions inside a repeated
            k-mer.
        examples: Up to five of the most frequent repeated k-mers, as
            ``(kmer, positions)``.
    """

    k: int
    redundant: int
    distinct_repeated: int
    max_occurrences: int
    covered_fraction: float
    examples: tuple[tuple[str, tuple[int, ...]], ...] = ()


def repeat_report(sequence: str, k: int = 12) -> RepeatReport:
    """Measure direct repeats of length ``k``.

    Unlike the original ``find_repeats``, which enumerated every substring
    length between ``min_length`` and ``max_length`` and so built a dictionary
    quadratic in sequence length, this indexes one k-mer size in a single
    linear pass. A repeat longer than ``k`` necessarily contains repeated
    k-mers, so it is still detected.

    Raises:
        ValueError: If ``k`` is not positive.
    """
    if k < 1:
        raise ValueError("k must be at least 1")
    positions: dict[str, list[int]] = {}
    for i in range(len(sequence) - k + 1):
        positions.setdefault(sequence[i : i + k], []).append(i)

    repeated = {kmer: pos for kmer, pos in positions.items() if len(pos) > 1}
    redundant = sum(len(pos) - 1 for pos in repeated.values())
    covered: set[int] = set()
    for kmer, pos in repeated.items():
        for start in pos:
            covered.update(range(start, start + k))
    examples = tuple(
        (kmer, tuple(pos))
        for kmer, pos in sorted(
            repeated.items(), key=lambda item: (-len(item[1]), item[0])
        )[:5]
    )
    return RepeatReport(
        k=k,
        redundant=redundant,
        distinct_repeated=len(repeated),
        max_occurrences=max((len(p) for p in positions.values()), default=0),
        covered_fraction=len(covered) / len(sequence) if sequence else 0.0,
        examples=examples,
    )


@dataclass(frozen=True)
class Hairpin:
    """An inverted repeat that can fold into a stem-loop.

    Attributes:
        stem_start: 0-based start of the 5' arm.
        stem_length: Length of the paired stem.
        loop_length: Number of unpaired bases between the arms.
        stem: The 5' arm sequence.
    """

    stem_start: int
    stem_length: int
    loop_length: int
    stem: str = field(default="", compare=False)

    @property
    def partner_start(self) -> int:
        """0-based start of the 3' arm."""
        return self.stem_start + self.stem_length + self.loop_length


def count_hairpin_pairs(
    sequence: str,
    *,
    stem: int = 10,
    min_loop: int = 3,
    max_loop: int = 60,
) -> int:
    """Count stem/loop arrangements with a stem of exactly ``stem`` bases.

    This is the quantity the optimizer penalises: one unit per pair of
    positions whose k-mers are reverse complements separated by an acceptable
    loop. A longer stem contains several such pairs and so is penalised more
    heavily, which is the desired behaviour.
    """
    _check_hairpin_params(stem, min_loop, max_loop)
    index: dict[str, list[int]] = {}
    for i in range(len(sequence) - stem + 1):
        index.setdefault(sequence[i : i + stem], []).append(i)

    pairs = 0
    for i in range(len(sequence) - stem + 1):
        partner = reverse_complement(sequence[i : i + stem])
        candidates = index.get(partner)
        if not candidates:
            continue
        low = i + stem + min_loop
        high = i + stem + max_loop
        pairs += bisect_right(candidates, high) - bisect_left(candidates, low)
    return pairs


def find_hairpins(
    sequence: str,
    *,
    stem: int = 10,
    min_loop: int = 3,
    max_loop: int = 60,
    limit: int = 50,
) -> list[Hairpin]:
    """Return hairpins with a stem of at least ``stem`` bases, stems extended.

    The original implementation looped over stem sizes, start positions and
    loop sizes with a substring search inside the innermost loop, which cost
    roughly cubic time (34 s for a 2.4 kb sequence in our benchmark). This
    indexes k-mers once and then extends each hit, which is linear in sequence
    length for a fixed loop range.
    """
    _check_hairpin_params(stem, min_loop, max_loop)
    if limit < 1:
        raise ValueError("limit must be positive")
    index: dict[str, list[int]] = {}
    for i in range(len(sequence) - stem + 1):
        index.setdefault(sequence[i : i + stem], []).append(i)

    seen: set[tuple[int, int]] = set()
    hairpins: list[Hairpin] = []
    for i in range(len(sequence) - stem + 1):
        partner = reverse_complement(sequence[i : i + stem])
        candidates = index.get(partner, ())
        first = bisect_left(candidates, i + stem + min_loop)
        last = bisect_right(candidates, i + stem + max_loop)
        for j in candidates[first:last]:  # j is the 3' arm start
            loop = j - (i + stem)
            if not (min_loop <= loop <= max_loop):
                continue
            left, right, length = i, j + stem - 1, stem
            # Extend the stem outwards while the flanking bases still pair.
            while (
                left > 0
                and right + 1 < len(sequence)
                and _pairs(sequence[left - 1], sequence[right + 1])
            ):
                left -= 1
                right += 1
                length += 1
            # Then inwards, so the same structure is described identically
            # however the k-mer scan happened to reach it.
            while (
                (right - left + 1) - 2 * length >= min_loop + 2
                and _pairs(sequence[left + length], sequence[right - length])
            ):
                length += 1
            key = (left, right)
            if key in seen:
                continue
            seen.add(key)
            hairpins.append(
                Hairpin(
                    stem_start=left,
                    stem_length=length,
                    loop_length=right - left + 1 - 2 * length,
                    stem=sequence[left : left + length],
                )
            )
            if len(hairpins) > limit:
                hairpins.sort(key=lambda h: (-h.stem_length, h.stem_start))
                del hairpins[limit:]
    hairpins.sort(key=lambda h: (-h.stem_length, h.stem_start))
    return hairpins


_PAIRS = frozenset({("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")})


def _pairs(left: str, right: str) -> bool:
    return (left, right) in _PAIRS


def _check_hairpin_params(stem: int, min_loop: int, max_loop: int) -> None:
    if stem < 1:
        raise ValueError("stem must be at least 1")
    if min_loop < 0 or max_loop < min_loop:
        raise ValueError("require 0 <= min_loop <= max_loop")


def motif_hits(sequence: str, motifs: list[Motif]) -> dict[str, tuple[int, ...]]:
    """Return the match positions of each motif that occurs in ``sequence``."""
    hits: dict[str, tuple[int, ...]] = {}
    for motif in motifs:
        found = motif.find_all(sequence)
        if found:
            hits[motif.name] = tuple(found)
    return hits
