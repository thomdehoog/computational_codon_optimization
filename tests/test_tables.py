"""Tests for codon usage tables and codon adaptation scoring."""

from __future__ import annotations

import json
import random

import pytest

from codonopt.sequence import AA_TO_CODONS
from codonopt.tables import (
    CodonUsage,
    available_tables,
    get_table,
    load_table,
)


@pytest.mark.parametrize("name", available_tables())
def test_bundled_tables_are_normalised(name):
    table = get_table(name)
    for amino_acid, codons in AA_TO_CODONS.items():
        total = sum(table.frequency(codon) for codon in codons)
        assert total == pytest.approx(1.0)
        assert max(table.weight(codon) for codon in codons) == pytest.approx(1.0)
        for codon in codons:
            assert 0.0 <= table.weight(codon) <= 1.0


def test_unknown_table_name_is_reported():
    with pytest.raises(KeyError, match="unknown codon usage table"):
        get_table("mouse")


def test_counts_are_accepted_and_normalised():
    table = CodonUsage.from_frequencies("counts", {"TTT": 30, "TTC": 70})
    assert table.frequency("TTT") == pytest.approx(0.3)
    assert table.weight("TTT") == pytest.approx(30 / 70)
    # Families absent from the input fall back to uniform rather than zero.
    assert table.frequency("GGT") == pytest.approx(0.25)


def test_missing_codons_within_a_family_are_never_chosen():
    table = CodonUsage.from_frequencies("partial", {"TTT": 1.0})
    rng = random.Random(0)
    assert {table.sample("F", rng) for _ in range(50)} == {"TTT"}


def test_from_frequencies_rejects_bad_input():
    with pytest.raises(ValueError, match="not a codon"):
        CodonUsage.from_frequencies("bad", {"XYZ": 1.0})
    with pytest.raises(ValueError, match="negative"):
        CodonUsage.from_frequencies("bad", {"TTT": -1.0})


def test_rna_spelling_is_accepted():
    table = CodonUsage.from_frequencies("rna", {"UUU": 1.0, "UUC": 1.0})
    assert table.frequency("TTT") == pytest.approx(0.5)


def test_cai_is_one_for_the_most_frequent_codons():
    table = get_table("human")
    best = "".join(table.most_frequent(aa) for aa in "MKWVTFISLL")
    assert table.cai(best) == pytest.approx(1.0)


def test_cai_penalises_rare_codons():
    table = get_table("human")
    protein = "LLLLLL"
    frequent = "".join(table.most_frequent(aa) for aa in protein)
    rare = "TTA" * len(protein)
    assert table.cai(rare) < table.cai(frequent)


def test_cai_ignores_single_codon_families():
    table = get_table("human")
    assert table.cai("ATGTGG") == 0.0


def test_sample_follows_the_distribution():
    table = CodonUsage.from_frequencies("skewed", {"AAA": 0.9, "AAG": 0.1})
    rng = random.Random(1)
    draws = [table.sample("K", rng) for _ in range(2000)]
    assert 0.85 < draws.count("AAA") / len(draws) < 0.95


def test_load_table_from_json_and_csv(tmp_path):
    payload = {"TTT": 1.0, "TTC": 3.0}
    json_path = tmp_path / "usage.json"
    json_path.write_text(json.dumps(payload))
    from_json = load_table(json_path)
    assert from_json.frequency("TTC") == pytest.approx(0.75)
    assert from_json.name == "usage"

    csv_path = tmp_path / "usage.csv"
    csv_path.write_text("codon,count\nTTT,1\nTTC,3\n")
    from_csv = load_table(csv_path, name="mine")
    assert from_csv.frequency("TTC") == pytest.approx(0.75)
    assert from_csv.name == "mine"


def test_load_table_rejects_unusable_files(tmp_path):
    bad_suffix = tmp_path / "usage.txt"
    bad_suffix.write_text("TTT 1")
    with pytest.raises(ValueError, match="expected a .json"):
        load_table(bad_suffix)

    empty = tmp_path / "empty.csv"
    empty.write_text("codon,count\n")
    with pytest.raises(ValueError, match="no codon rows"):
        load_table(empty)

    malformed = tmp_path / "malformed.csv"
    malformed.write_text("TTT,abc\n")
    with pytest.raises(ValueError, match="not a number"):
        load_table(malformed)


def test_table_from_text_matches_load_table(tmp_path):
    from codonopt.tables import table_from_text

    text = "codon,count\nTTT,1\nTTC,3\n"
    path = tmp_path / "usage.csv"
    path.write_text(text)
    from_text = table_from_text(text, "csv", name="usage")
    from_file = load_table(path)
    assert from_text.frequencies == from_file.frequencies
    assert table_from_text("TTT\t1\nTTC\t3\n", ".tsv", name="t").frequency("TTC") == 0.75
    assert table_from_text('{"TTT": 1, "TTC": 3}', "json", name="j").frequency("TTC") == 0.75


def test_table_from_text_errors_name_their_source():
    from codonopt.tables import table_from_text

    with pytest.raises(ValueError, match="pasted:2: 'x' is not a number"):
        table_from_text("TTT,1\nTTC,x\n", "csv", name="n", source="pasted")
    with pytest.raises(ValueError, match="pasted: every JSON value must be a number"):
        table_from_text('{"TTT": "high"}', "json", name="n", source="pasted")
    with pytest.raises(ValueError, match="pasted: expected a JSON object"):
        table_from_text("[1, 2]", "json", name="n", source="pasted")
    with pytest.raises(ValueError, match="expected a .json, .csv or .tsv"):
        table_from_text("TTT 1", "txt", name="n", source="pasted")
