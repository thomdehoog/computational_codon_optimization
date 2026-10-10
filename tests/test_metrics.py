"""Tests for the measurement functions."""

from __future__ import annotations

import random

import pytest

from codonopt.metrics import (
    count_hairpin_pairs,
    find_hairpins,
    gc_window_report,
    homopolymer_excess,
    homopolymer_runs,
    motif_hits,
    repeat_report,
    window_starts,
)
from codonopt.motifs import resolve_motifs
from codonopt.sequence import reverse_complement


def test_window_starts_cover_the_tail():
    starts = window_starts(100, 30, 15)
    assert starts[0] == 0
    assert starts[-1] == 70  # the final window ends exactly at the sequence end
    assert all(start + 30 <= 100 for start in starts)
    assert window_starts(10, 30, 3) == [0]  # sequence shorter than one window
    assert window_starts(0, 30, 3) == []


def test_window_starts_validates_arguments():
    with pytest.raises(ValueError):
        window_starts(100, 0, 3)


def test_gc_window_report_counts_both_directions():
    report = gc_window_report(
        "A" * 30 + "G" * 30, window=30, step=30, target_range=(0.4, 0.6)
    )
    assert report.low == 1
    assert report.high == 1
    assert report.violations == 2
    assert report.overall == pytest.approx(0.5)
    assert report.min_fraction == pytest.approx(0.0)
    assert report.max_fraction == pytest.approx(1.0)


def test_gc_window_report_rejects_percentages():
    with pytest.raises(ValueError, match="use 0.3 rather than 30"):
        gc_window_report("ACGT" * 10, target_range=(30, 70))


def test_gc_window_report_handles_empty_input():
    report = gc_window_report("")
    assert report.windows == 0
    assert report.violations == 0


def test_homopolymer_runs_and_excess():
    sequence = "AC" + "T" * 9 + "GG"
    runs = homopolymer_runs(sequence, 3)
    assert [(run.base, run.start, run.length) for run in runs] == [("T", 2, 9)]
    # A run of 9 with a limit of 6 pushes three bases past the limit.
    assert homopolymer_excess(sequence, 6) == 3
    assert homopolymer_excess(sequence, 9) == 0


def test_homopolymer_excess_matches_run_definition():
    rng = random.Random(4)
    sequence = "".join(rng.choice("AACCGGTT") for _ in range(400))
    for limit in (3, 4, 6):
        expected = sum(
            run.length - limit for run in homopolymer_runs(sequence, limit + 1)
        )
        assert homopolymer_excess(sequence, limit) == expected


def test_repeat_report_finds_a_planted_repeat():
    unit = "ACGTTGCAACGT"
    sequence = unit + "AAACCCGGGTTT" + unit
    report = repeat_report(sequence, k=len(unit))
    assert report.redundant == 1
    assert report.distinct_repeated == 1
    assert report.max_occurrences == 2
    assert report.examples[0][0] == unit


def test_repeat_report_is_zero_for_unique_sequence():
    rng = random.Random(5)
    sequence = "".join(rng.choice("ACGT") for _ in range(600))
    assert repeat_report(sequence, k=14).redundant == 0


def test_repeat_report_validates_k():
    with pytest.raises(ValueError):
        repeat_report("ACGT", k=0)


def test_hairpin_detection_finds_a_planted_stem_loop():
    stem = "ACGTTGCACC"
    loop = "TTTTT"
    sequence = "GGCAT" + stem + loop + reverse_complement(stem) + "ATTAC"
    hairpins = find_hairpins(sequence, stem=len(stem), min_loop=3, max_loop=20)
    assert hairpins
    best = hairpins[0]
    assert best.stem_length >= len(stem)
    assert best.loop_length == len(loop)
    # The flanks extend the stem to 12 bp, which holds three overlapping
    # 10-mer pairs; a longer stem is meant to be penalised more heavily.
    assert best.stem_length == 12
    assert count_hairpin_pairs(sequence, stem=len(stem), min_loop=3, max_loop=20) == 3


def test_hairpin_pairs_respects_the_loop_range():
    stem = "ACGTTGCACC"
    sequence = stem + "T" * 40 + reverse_complement(stem)
    assert count_hairpin_pairs(sequence, stem=10, min_loop=3, max_loop=20) == 0
    assert count_hairpin_pairs(sequence, stem=10, min_loop=3, max_loop=60) == 1


def test_hairpin_arguments_are_validated():
    with pytest.raises(ValueError):
        count_hairpin_pairs("ACGT", stem=4, min_loop=5, max_loop=2)


def test_motif_hits_reports_both_strands():
    motifs = resolve_motifs(["EcoRI"])
    assert motif_hits("TTGAATTCTT", motifs) == {"EcoRI": (2,)}
    # The reverse complement of GAATTC is itself, so a one-sided site needs
    # an asymmetric motif to test.
    asymmetric = resolve_motifs(["BsaI"])  # GGTCTC / GAGACC
    assert motif_hits("AAGAGACCAA", asymmetric) == {"BsaI": (2,)}
    assert motif_hits("AAAAAA", asymmetric) == {}


def test_hairpin_description_is_register_independent():
    """The same structure must be reported once, with its maximal stem.

    A stem/loop can be reached from several k-mer registers. Without
    extending the stem inwards as well as outwards, whichever register the
    scan reached first decided the reported stem length and loop size.
    """
    stem = "ACGTTGCACC"
    sequence = "GGCAT" + stem + "TTTTT" + reverse_complement(stem) + "ATTAC"
    for k in (8, 9, 10):
        hairpins = find_hairpins(sequence, stem=k, min_loop=3, max_loop=20)
        best = hairpins[0]
        assert (best.stem_start, best.stem_length, best.loop_length) == (3, 12, 5)


def test_find_hairpins_respects_the_limit():
    unit = "ACGTTGCACC"
    sequence = (unit + "TTTTT" + reverse_complement(unit) + "AAA") * 6
    assert len(find_hairpins(sequence, stem=10, min_loop=3, max_loop=20, limit=3)) == 3
