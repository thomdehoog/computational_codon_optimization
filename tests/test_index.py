"""The incremental k-mer index must agree with a full recomputation.

These are the tests that guard the optimizer's inner loop: if the index drifts
from the brute-force definitions in :mod:`codonopt.metrics`, every score the
optimizer computes is wrong, and nothing else would notice.
"""

from __future__ import annotations

import random

import pytest

from codonopt.index import KmerIndex
from codonopt.metrics import count_hairpin_pairs, repeat_report
from codonopt.sequence import reverse_complement

K = 8
MIN_LOOP = 3
MAX_LOOP = 24


def build(sequence: str) -> KmerIndex:
    return KmerIndex(
        sequence, K, min_loop=MIN_LOOP, max_loop=MAX_LOOP, track_hairpins=True
    )


def brute_force(sequence: str) -> tuple[int, int]:
    return (
        repeat_report(sequence, K).redundant,
        count_hairpin_pairs(sequence, stem=K, min_loop=MIN_LOOP, max_loop=MAX_LOOP),
    )


def test_fresh_index_matches_brute_force():
    rng = random.Random(0)
    sequence = "".join(rng.choice("ACGT") for _ in range(300))
    index = build(sequence)
    expected = brute_force(sequence)
    assert (index.redundant, index.hairpin_pairs) == expected


def test_index_counts_planted_features():
    unit = "ACGTTGCA"
    sequence = unit + "TTTTT" + unit  # a direct repeat
    index = build(sequence)
    assert index.redundant == 1
    assert index.occurrences(unit) == 2

    # One stem/loop arrangement is one pair, counted once rather than once
    # per endpoint.
    hairpin = unit + "TTTTT" + reverse_complement(unit)
    assert build(hairpin).hairpin_pairs == 1


def test_incremental_updates_stay_exact_under_random_edits():
    rng = random.Random(1)
    # A sequence built from a small alphabet of blocks, so repeats and
    # hairpins actually appear and then disappear as edits are applied.
    blocks = ["ACGT", "ACGA", "TTGC", "GCAA", "TGCA"]
    sequence = "".join(rng.choice(blocks) for _ in range(60))
    index = build(sequence)

    for _ in range(200):
        position = rng.randrange(0, len(sequence) - 3)
        replacement = "".join(rng.choice("ACGT") for _ in range(3))
        sequence = sequence[:position] + replacement + sequence[position + 3 :]
        index.replace_range(sequence, position, position + 3)
        expected = brute_force(sequence)
        assert (index.redundant, index.hairpin_pairs) == expected


def test_incremental_updates_at_the_sequence_edges():
    rng = random.Random(2)
    sequence = "".join(rng.choice("ACGT") for _ in range(60))
    index = build(sequence)
    for position in (0, 3, len(sequence) - 3, len(sequence) - 6):
        replacement = "".join(rng.choice("ACGT") for _ in range(3))
        sequence = sequence[:position] + replacement + sequence[position + 3 :]
        index.replace_range(sequence, position, position + 3)
        assert (index.redundant, index.hairpin_pairs) == brute_force(sequence)


def test_remove_then_add_restores_the_counters():
    rng = random.Random(3)
    sequence = "".join(rng.choice("ACG") for _ in range(80))
    index = build(sequence)
    before = (index.redundant, index.hairpin_pairs, index.distinct)
    index.replace_range(sequence, 30, 33)  # same bases back again
    assert (index.redundant, index.hairpin_pairs, index.distinct) == before


def test_affected_positions_cover_every_overlapping_kmer():
    index = KmerIndex("A" * 20, 5)
    assert index.affected_positions(10, 13) == range(6, 13)
    assert index.affected_positions(0, 3) == range(0, 3)
    # Near the end, positions past the last k-mer start are dropped.
    assert index.affected_positions(17, 20) == range(13, 16)


def test_sequence_shorter_than_k_is_empty():
    index = KmerIndex("ACGT", 8)
    assert index.span == 0
    assert index.redundant == 0
    assert index.hairpin_pairs == 0
    assert index.kmer_at(0) is None


def test_invalid_parameters_are_rejected():
    with pytest.raises(ValueError):
        KmerIndex("ACGT", 0)
    with pytest.raises(ValueError):
        KmerIndex("ACGT", 2, min_loop=5, max_loop=1)


def test_hairpin_tracking_can_be_switched_off():
    unit = "ACGTTGCA"
    sequence = unit + "TTTTT" + reverse_complement(unit)
    index = KmerIndex(sequence, len(unit), track_hairpins=False)
    assert index.hairpin_pairs == 0
    assert index.redundant == 0
