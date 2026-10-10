"""An incrementally maintained k-mer index.

Direct repeats and hairpins are the two global properties of a coding sequence
that the optimizer has to watch, and both are expensive to recompute from
scratch after every candidate codon swap. This index keeps the two quantities
up to date under local edits:

``redundant``
    Occurrences of a k-mer beyond its first, summed over all k-mers. Zero means
    every k-mer in the sequence is unique.

``hairpin_pairs``
    Pairs of positions whose k-mers are reverse complements of each other,
    separated by a loop of an acceptable length. Each pair is counted once,
    which matches :func:`codonopt.metrics.count_hairpin_pairs`. A stem longer
    than ``k`` contains several overlapping pairs and so counts for more.

Removing a position subtracts exactly the pairs that position takes part in,
and adding one adds the pairs it forms against the positions already present,
so the counters stay exact however many edits are applied.
"""

from __future__ import annotations

from .sequence import reverse_complement

__all__ = ["KmerIndex"]


class KmerIndex:
    """Index of every k-mer start position in a sequence.

    Args:
        sequence: The sequence to index.
        k: K-mer size.
        min_loop: Smallest loop accepted between hairpin arms.
        max_loop: Largest loop accepted between hairpin arms.
        track_hairpins: Set False to skip pair counting when the hairpin term
            is switched off, which makes updates noticeably cheaper.

    Raises:
        ValueError: On a non-positive ``k`` or an invalid loop range.
    """

    __slots__ = (
        "k",
        "min_loop",
        "max_loop",
        "track_hairpins",
        "_kmer_at",
        "_counts",
        "_total",
        "_distinct",
        "_pairs",
        "_offsets",
    )

    def __init__(
        self,
        sequence: str,
        k: int,
        *,
        min_loop: int = 3,
        max_loop: int = 60,
        track_hairpins: bool = True,
    ) -> None:
        if k < 1:
            raise ValueError("k must be at least 1")
        if min_loop < 0 or max_loop < min_loop:
            raise ValueError("require 0 <= min_loop <= max_loop")
        self.k = k
        self.min_loop = min_loop
        self.max_loop = max_loop
        self.track_hairpins = track_hairpins
        self._offsets = range(k + min_loop, k + max_loop + 1)
        self._counts: dict[str, int] = {}
        self._total = 0
        self._distinct = 0
        self._pairs = 0
        slots = max(0, len(sequence) - k + 1)
        self._kmer_at: list[str | None] = [None] * slots
        for position in range(slots):
            self.add(position, sequence[position : position + k])

    @property
    def redundant(self) -> int:
        """Number of k-mer occurrences beyond the first of each distinct k-mer."""
        return self._total - self._distinct

    @property
    def distinct(self) -> int:
        """Number of distinct k-mers currently indexed."""
        return self._distinct

    @property
    def hairpin_pairs(self) -> int:
        """Number of reverse-complement pairs within the loop range.

        Each pair of positions is counted once, not once per endpoint.
        """
        return self._pairs

    @property
    def span(self) -> int:
        """Number of k-mer start positions in the indexed sequence."""
        return len(self._kmer_at)

    def kmer_at(self, position: int) -> str | None:
        """Return the k-mer starting at ``position``, if one is indexed."""
        if 0 <= position < len(self._kmer_at):
            return self._kmer_at[position]
        return None

    def occurrences(self, kmer: str) -> int:
        """Return how many times ``kmer`` occurs."""
        return self._counts.get(kmer, 0)

    def affected_positions(self, lo: int, hi: int) -> range:
        """Return the k-mer start positions that overlap bases ``[lo, hi)``."""
        start = max(0, lo - self.k + 1)
        stop = min(len(self._kmer_at), hi)
        return range(start, max(start, stop))

    def add(self, position: int, kmer: str) -> None:
        """Index ``kmer`` at ``position``."""
        self._kmer_at[position] = kmer
        count = self._counts.get(kmer, 0)
        if count == 0:
            self._distinct += 1
        self._counts[kmer] = count + 1
        self._total += 1
        if self.track_hairpins:
            self._pairs += self._pairs_at(position, kmer)

    def remove(self, position: int) -> None:
        """Drop whatever k-mer is indexed at ``position``."""
        kmer = self._kmer_at[position]
        if kmer is None:
            return
        if self.track_hairpins:
            self._pairs -= self._pairs_at(position, kmer)
        count = self._counts[kmer] - 1
        if count:
            self._counts[kmer] = count
        else:
            del self._counts[kmer]
            self._distinct -= 1
        self._total -= 1
        self._kmer_at[position] = None

    def _pairs_at(self, position: int, kmer: str) -> int:
        """Count hairpin pairs that ``position`` takes part in right now."""
        partner = reverse_complement(kmer)
        if partner not in self._counts:
            return 0
        kmer_at = self._kmer_at
        limit = len(kmer_at)
        pairs = 0
        for offset in self._offsets:
            downstream = position + offset
            if downstream < limit and kmer_at[downstream] == partner:
                pairs += 1
            upstream = position - offset
            if upstream >= 0 and kmer_at[upstream] == partner:
                pairs += 1
        return pairs

    def replace_range(self, sequence: str, lo: int, hi: int) -> None:
        """Re-index the positions overlapping ``[lo, hi)`` from ``sequence``.

        ``sequence`` must already contain the new bases.
        """
        positions = self.affected_positions(lo, hi)
        for position in positions:
            self.remove(position)
        k = self.k
        for position in positions:
            self.add(position, sequence[position : position + k])
