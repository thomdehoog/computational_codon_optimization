"""The incrementally maintained score must equal a full recomputation."""

from __future__ import annotations

import random

import pytest

from codonopt import DesignConfig, Weights, get_table
from codonopt.motifs import resolve_motifs
from codonopt.sequence import CODON_TO_AA, translate
from codonopt.splice import MotifScorer
from codonopt.state import SequenceState

from conftest import duplicated_pair, duplicated_protein, random_dna, random_protein


def test_score_matches_full_recomputation_after_random_swaps(config):
    protein = duplicated_protein(block=40, seed=1)
    state = SequenceState(random_dna(protein, seed=2), config)
    rng = random.Random(3)

    for _ in range(120):
        position = rng.randrange(state.codon_count)
        options = state.alternatives(position)
        state.set_codon(position, rng.choice(options))
        assert state.score == pytest.approx(state.audit(), abs=1e-9)


def test_breakdown_total_equals_the_running_score(config):
    protein = random_protein(80, seed=4)
    state = SequenceState(random_dna(protein, seed=5), config)
    rng = random.Random(6)
    for _ in range(25):
        position = rng.randrange(state.codon_count)
        state.set_codon(position, rng.choice(state.alternatives(position)))
    assert state.breakdown().total == pytest.approx(state.score, abs=1e-9)


def test_every_term_is_exercised_by_the_test_configuration(config):
    _, dna = duplicated_pair(block=40, seed=7)
    state = SequenceState(dna, config)
    breakdown = state.breakdown()
    assert set(breakdown.terms) == {
        "gc_global",
        "gc_window",
        "repeat",
        "hairpin",
        "homopolymer",
        "motif",
        "splice",
        "usage",
    }
    # The duplicated coding sequence guarantees the repeat term is non-zero.
    assert breakdown.raw["repeat"] > 0


def test_setting_the_same_codon_is_a_no_op(config):
    state = SequenceState(random_dna(random_protein(30, seed=9), seed=10), config)
    before = state.score
    sequence = state.sequence
    assert state.set_codon(3, state.codon(3)) == state.codon(3)
    assert state.sequence == sequence
    assert state.score == before


def test_non_synonymous_replacement_is_refused(config):
    state = SequenceState("ATGAAATGG", config)
    with pytest.raises(ValueError, match="encodes"):
        state.set_codon(1, "GGG")
    assert translate(state.sequence) == "MKW"


def test_set_codon_preserves_the_protein(config):
    protein = random_protein(60, seed=11)
    state = SequenceState(random_dna(protein, seed=12), config)
    rng = random.Random(13)
    for _ in range(100):
        position = rng.randrange(state.codon_count)
        state.set_codon(position, rng.choice(state.alternatives(position)))
    assert translate(state.sequence) == protein
    assert state.protein == protein


def test_motif_penalty_counts_literal_and_iupac_patterns():
    # GGTCTC (BsaI) on the forward strand and an IUPAC pattern that matches it.
    config = DesignConfig(
        table=get_table("uniform"),
        splice_scorer=None,
        motifs=tuple(resolve_motifs(["BsaI", "ambig=GGTCTN"])),
        weights=Weights(
            gc_global=0, gc_window=0, repeat=0, hairpin=0, homopolymer=0,
            motif=1.0, splice=0, usage=0,
        ),
    )
    with_site = SequenceState("ATGGGTCTCAAA", config)
    # One literal hit plus one IUPAC hit at the same position.
    assert with_site.score == pytest.approx(2.0)
    assert SequenceState("ATGAAAAAACCC", config).score == pytest.approx(0.0)


def test_homopolymer_term_counts_bases_past_the_limit():
    config = DesignConfig(
        table=get_table("uniform"),
        splice_scorer=None,
        max_homopolymer=4,
        weights=Weights(
            gc_global=0, gc_window=0, repeat=0, hairpin=0, homopolymer=1.0,
            motif=0, splice=0, usage=0,
        ),
    )
    # A run of eight with a limit of four leaves four bases over the limit.
    state = SequenceState("AAAAAAAA" + "C" + "ATG", config)
    assert state.score == pytest.approx(4.0)


def test_usage_term_prefers_frequent_codons():
    table = get_table("human")
    config = DesignConfig(
        table=table,
        splice_scorer=None,
        weights=Weights(
            gc_global=0, gc_window=0, repeat=0, hairpin=0, homopolymer=0,
            motif=0, splice=0, usage=1.0,
        ),
    )
    frequent = SequenceState("ATG" + table.most_frequent("L") * 5, config)
    rare = SequenceState("ATG" + min(
        ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"), key=table.weight
    ) * 5, config)
    assert frequent.score < rare.score


def test_splice_term_responds_to_a_strong_donor():
    config = DesignConfig(
        table=get_table("uniform"),
        splice_scorer=MotifScorer(),
        weights=Weights(
            gc_global=0, gc_window=0, repeat=0, hairpin=0, homopolymer=0,
            motif=0, splice=1.0, usage=0,
        ),
    )
    # CAG|GTAAGT is the canonical donor the motif screen scores highest.
    strong = SequenceState("CAGGTAAGTATG", config)
    weak = SequenceState("ATGAAACCCTTT", config)
    assert strong.score > weak.score


def test_score_is_stable_across_long_edit_chains(config):
    """A long run of edits must not accumulate floating point drift."""
    protein = random_protein(100, seed=14)
    state = SequenceState(random_dna(protein, seed=15), config)
    rng = random.Random(16)
    for _ in range(500):
        position = rng.randrange(state.codon_count)
        state.set_codon(position, rng.choice(state.alternatives(position)))
    assert state.score == pytest.approx(state.audit(), abs=1e-7)


def test_alternatives_are_synonymous(config):
    state = SequenceState(random_dna(random_protein(40, seed=17), seed=18), config)
    for position in range(state.codon_count):
        residue = CODON_TO_AA[state.codon(position)]
        assert all(
            CODON_TO_AA[option] == residue for option in state.alternatives(position)
        )


def test_breakdown_agrees_with_the_score_for_overlapping_motifs():
    """A site must cost the same in the objective as it does in the report.

    Patterns that are reverse complements of each other, and ambiguous
    patterns that match both strands at one position, are the cases where a
    per-position charge and a per-motif count can disagree.
    """
    only_motifs = Weights(
        gc_global=0, gc_window=0, repeat=0, hairpin=0, homopolymer=0,
        motif=1.0, splice=0, usage=0,
    )
    for specs, sequence, expected in (
        (["GAAGAC", "GTCTTC"], "ATGGAAGACTGG", 1.0),
        (["x=NGATCY"], "GGATCCATG", 1.0),
        (["BsaI", "ambig=GGTCTN"], "ATGGGTCTCAAA", 2.0),
        (["EcoRI"], "ATGGAATTCAAA", 1.0),
    ):
        config = DesignConfig(
            table=get_table("uniform"),
            splice_scorer=None,
            motifs=tuple(resolve_motifs(specs)),
            weights=only_motifs,
        )
        state = SequenceState(sequence, config)
        assert state.score == pytest.approx(expected), specs
        assert state.breakdown().total == pytest.approx(state.score), specs
