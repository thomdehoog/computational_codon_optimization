"""Tests for motif resolution and matching."""

from __future__ import annotations

import pytest

from codonopt.motifs import (
    MOTIF_GROUPS,
    REGULATORY_MOTIFS,
    RESTRICTION_SITES,
    Motif,
    iupac_to_regex,
    resolve_motifs,
)
from codonopt.sequence import reverse_complement


def test_every_bundled_pattern_is_valid_iupac():
    for table in (RESTRICTION_SITES, REGULATORY_MOTIFS):
        for name, pattern in table.items():
            assert pattern == pattern.upper()
            iupac_to_regex(pattern)  # raises if a code is unknown
            assert Motif(name, pattern).length == len(pattern)


def test_groups_only_reference_known_names():
    for group, members in MOTIF_GROUPS.items():
        for member in members:
            assert member in RESTRICTION_SITES or member in REGULATORY_MOTIFS


def test_iupac_codes_expand_to_character_classes():
    assert iupac_to_regex("ACGT") == "ACGT"
    assert iupac_to_regex("RYN") == "[AG][CT][ACGT]"
    with pytest.raises(ValueError, match="not an IUPAC"):
        iupac_to_regex("ACZ")


def test_motif_matches_both_strands():
    motif = Motif("BsaI", "GGTCTC")
    assert motif.find_all("AAGGTCTCAA") == [2]
    assert motif.find_all("AA" + reverse_complement("GGTCTC") + "AA") == [2]
    assert motif.find_all("AAAAAA") == []


def test_palindromic_motif_is_counted_once():
    motif = Motif("EcoRI", "GAATTC")
    assert motif.find_all("AAGAATTCAA") == [2]


def test_overlapping_matches_are_all_reported():
    motif = Motif("poly", "AATAAA")
    # AATAAATAAA holds two overlapping copies.
    assert motif.find_all("AATAAATAAA") == [0, 4]


def test_iupac_motif_matches_every_variant():
    motif = Motif("ambiguous", "GGTCTN")
    for base in "ACGT":
        assert motif.find_all("AA" + "GGTCT" + base + "AA") == [2]


def test_resolve_accepts_groups_names_patterns_and_assignments():
    motifs = resolve_motifs(["golden_gate", "EcoRI", "AATAAA", "mine=GGWCC"])
    patterns = {motif.pattern for motif in motifs}
    assert RESTRICTION_SITES["BsaI"] in patterns
    assert RESTRICTION_SITES["EcoRI"] in patterns
    assert "AATAAA" in patterns
    assert "GGWCC" in patterns
    assert any(motif.name == "mine" for motif in motifs)


def test_resolve_is_case_insensitive_for_names():
    assert resolve_motifs(["ecori"])[0].pattern == RESTRICTION_SITES["EcoRI"]
    assert resolve_motifs(["ecori"])[0].name == "EcoRI"


def test_resolve_collapses_duplicates():
    motifs = resolve_motifs(["EcoRI", "GAATTC", "common_cloning"])
    patterns = [motif.pattern for motif in motifs]
    assert len(patterns) == len(set(patterns))


def test_resolve_handles_empty_input():
    assert resolve_motifs(None) == []
    assert resolve_motifs([]) == []
    assert resolve_motifs(["  "]) == []
    assert resolve_motifs(["none"]) == []


def test_resolve_rejects_nonsense():
    with pytest.raises(ValueError, match="cannot resolve motif"):
        resolve_motifs(["NotAnEnzyme"])


def test_empty_pattern_is_refused():
    with pytest.raises(ValueError, match="empty pattern"):
        Motif("blank", "   ")


def test_reverse_complement_patterns_name_one_site():
    """Both strands are searched anyway, so the pair must not count twice."""
    motifs = resolve_motifs(["GAAGAC", "GTCTTC"])
    assert len(motifs) == 1
    assert motifs[0].find_all("ATGGAAGACTGG") == [3]
    # Resolution order must not change the outcome.
    assert len(resolve_motifs(["GTCTTC", "GAAGAC"])) == 1
    # BbsI is GAAGAC, so naming the enzyme and its other strand is also one.
    assert len(resolve_motifs(["BbsI", "GTCTTC"])) == 1


def test_a_position_matching_both_strands_is_one_site():
    # NGATCY matches GGATCC, and so does its reverse complement RGATCN.
    motif = Motif("ambiguous", "NGATCY")
    assert motif.find_all("GGATCCATG") == [0]
