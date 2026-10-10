"""Tests for the nucleotide and protein primitives."""

from __future__ import annotations

import random

import pytest

from codonopt.sequence import (
    AA_TO_CODONS,
    CODON_TO_AA,
    SequenceError,
    gc_fraction,
    gc_prefix,
    iter_codons,
    reverse_complement,
    split_codons,
    translate,
    validate_dna,
    validate_protein,
    window_gc_count,
)


def test_genetic_code_is_complete_and_consistent():
    assert len(CODON_TO_AA) == 64
    for amino_acid, codons in AA_TO_CODONS.items():
        assert len(set(codons)) == len(codons)
        for codon in codons:
            assert len(codon) == 3
            assert CODON_TO_AA[codon] == amino_acid


def test_translate_round_trips_every_codon():
    for codon, amino_acid in CODON_TO_AA.items():
        assert translate(codon) == amino_acid


def test_translate_rejects_partial_codons_and_junk():
    with pytest.raises(SequenceError, match="multiple of three"):
        translate("ATGA")
    with pytest.raises(SequenceError, match="not a valid codon"):
        translate("ATGNNN")


def test_reverse_complement_is_an_involution():
    rng = random.Random(0)
    sequence = "".join(rng.choice("ACGT") for _ in range(99))
    assert reverse_complement(reverse_complement(sequence)) == sequence
    assert reverse_complement("AACGT") == "ACGTT"


def test_gc_helpers_agree():
    sequence = "GGCCATATAT"
    assert gc_fraction(sequence) == pytest.approx(0.4)
    prefix = gc_prefix(sequence)
    assert window_gc_count(prefix, 0, 4) == 4
    assert window_gc_count(prefix, 4, len(sequence)) == 0
    assert gc_fraction("") == 0.0


def test_split_and_iter_codons():
    assert split_codons("ATGAAA") == ["ATG", "AAA"]
    assert list(iter_codons("ATGAAAG")) == ["ATG", "AAA"]
    with pytest.raises(SequenceError):
        split_codons("ATGAA")


@pytest.mark.parametrize(
    "sequence, message",
    [("", "empty"), ("ATGX", "non-ACGT"), ("atgaa", "multiple of three")],
)
def test_validate_dna_rejects_bad_input(sequence, message):
    with pytest.raises(SequenceError, match=message):
        validate_dna(sequence, in_frame=True)


def test_validate_dna_normalises_case_and_whitespace():
    assert validate_dna("  atgaaa\n") == "ATGAAA"


def test_validate_protein_allows_trailing_stop_only():
    assert validate_protein("MKV*") == "MKV*"
    with pytest.raises(SequenceError, match="internal stop"):
        validate_protein("MK*V")
    assert validate_protein("MK*V", allow_internal_stop=True) == "MK*V"
    with pytest.raises(SequenceError, match="unsupported residues"):
        validate_protein("MKZ")
