"""Tests for FASTA reading and writing."""

from __future__ import annotations

import pytest

from codonopt.fasta import Record, format_fasta, parse_fasta, read_fasta, write_fasta


def test_parse_single_record():
    records = list(parse_fasta([">gene one", "ATGAAA", "GGG"]))
    assert len(records) == 1
    assert records[0].name == "gene"
    assert records[0].description == "gene one"
    assert records[0].sequence == "ATGAAAGGG"


def test_parse_multiple_records_and_skip_comments():
    text = [">a", "ATG", "; comment", "# comment", ">b", "AAA", "", "CCC"]
    records = list(parse_fasta(text))
    assert [record.name for record in records] == ["a", "b"]
    assert records[1].sequence == "AAACCC"


def test_parse_bare_sequence_without_header():
    records = list(parse_fasta(["ATGAAA"]))
    assert records[0].name == "sequence"
    assert records[0].sequence == "ATGAAA"


def test_parse_strips_internal_spaces():
    records = list(parse_fasta([">a", "ATG AAA"]))
    assert records[0].sequence == "ATGAAA"


def test_header_without_sequence_is_an_error():
    with pytest.raises(ValueError, match="no sequence"):
        list(parse_fasta([">lonely"]))


def test_format_wraps_at_the_requested_width():
    record = Record("x", "x description", "A" * 25)
    text = format_fasta(record, width=10)
    lines = text.splitlines()
    assert lines[0] == ">x description"
    assert lines[1:] == ["A" * 10, "A" * 10, "A" * 5]
    with pytest.raises(ValueError):
        format_fasta(record, width=0)


def test_round_trip_through_a_file(tmp_path):
    path = tmp_path / "out.fna"
    records = [Record("a", "a first", "ATGAAA"), Record("b", "b second", "GGGCCC")]
    write_fasta(path, records, width=4)
    loaded = read_fasta(path)
    assert [record.sequence for record in loaded] == ["ATGAAA", "GGGCCC"]
    assert [record.name for record in loaded] == ["a", "b"]


def test_read_empty_file_is_an_error(tmp_path):
    path = tmp_path / "empty.fna"
    path.write_text("\n\n")
    with pytest.raises(ValueError, match="no FASTA records"):
        read_fasta(path)
