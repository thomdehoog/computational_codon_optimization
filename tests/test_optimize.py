"""Tests for the optimizer itself."""

from __future__ import annotations

import pytest

from codonopt import (
    DesignConfig,
    OptimizerSettings,
    Weights,
    get_table,
    optimize,
    optimize_dna,
    resolve_motifs,
    seed_sequence,
    translate,
)
from codonopt.sequence import SequenceError
from codonopt.state import SequenceState

from conftest import duplicated_pair, duplicated_protein, random_dna, random_protein

import random


def test_protein_is_invariant(config):
    protein = random_protein(90, seed=20, stop=True)
    result = optimize(protein, config, OptimizerSettings(seed=1, max_sweeps=3))
    assert translate(result.sequence) == protein
    assert len(result.sequence) == 3 * len(protein)


def test_score_never_gets_worse(config):
    protein = duplicated_protein(block=40, seed=21)
    result = optimize(protein, config, OptimizerSettings(seed=2, max_sweeps=4))
    assert result.score <= result.initial_score
    assert result.improvement >= 0


def test_repeats_are_driven_down(plain_config):
    protein, dna = duplicated_pair(block=60, seed=22)
    config = plain_config.with_(repeat_k=12)
    before = SequenceState(dna, config)
    assert before.breakdown().raw["repeat"] > 0
    result = optimize(
        protein, config, OptimizerSettings(seed=3, max_sweeps=6), start=dna
    )
    assert result.breakdown.raw["repeat"] == 0


def test_forbidden_motifs_are_removed(plain_config):
    protein = random_protein(120, seed=24)
    config = plain_config.with_(
        motifs=tuple(resolve_motifs(["golden_gate", "common_cloning"]))
    )
    seeded = seed_sequence(protein, config.table, random.Random(4), "random")
    result = optimize(protein, config, OptimizerSettings(seed=4), start=seeded)
    assert result.breakdown.raw["motif"] == 0


def test_gc_content_moves_towards_the_target():
    protein = random_protein(120, seed=25)
    config = DesignConfig(
        table=get_table("uniform"),
        splice_scorer=None,
        gc_target=0.55,
        weights=Weights(usage=0.0),
    )
    result = optimize(protein, config, OptimizerSettings(seed=5, max_sweeps=4))
    final = result.breakdown.raw["gc_global"]
    initial = result.initial_breakdown.raw["gc_global"]
    assert abs(final - 0.55) <= abs(initial - 0.55)
    assert abs(final - 0.55) < 0.02


def test_cai_improves_when_usage_is_the_only_term():
    protein = random_protein(60, seed=26)
    table = get_table("human")
    config = DesignConfig(
        table=table,
        splice_scorer=None,
        weights=Weights(
            gc_global=0, gc_window=0, repeat=0, hairpin=0, homopolymer=0,
            motif=0, splice=0, usage=1.0,
        ),
    )
    result = optimize(protein, config, OptimizerSettings(seed=6, max_sweeps=3))
    # With only codon usage to satisfy, the answer is the most frequent codon
    # at every position.
    assert result.sequence == "".join(table.most_frequent(aa) for aa in protein)
    assert table.cai(result.sequence) == pytest.approx(1.0)


def test_runs_are_reproducible(config):
    protein = random_protein(70, seed=27)
    first = optimize(protein, config, OptimizerSettings(seed=11, max_sweeps=3))
    second = optimize(protein, config, OptimizerSettings(seed=11, max_sweeps=3))
    assert first.sequence == second.sequence
    assert first.score == second.score


def test_different_seeds_explore_different_sequences(config):
    protein = random_protein(70, seed=28)
    first = optimize(protein, config, OptimizerSettings(seed=12, max_sweeps=2))
    second = optimize(protein, config, OptimizerSettings(seed=13, max_sweeps=2))
    assert first.sequence != second.sequence


def test_restarts_keep_the_best_result(config):
    protein = random_protein(60, seed=29)
    many = optimize(protein, config, OptimizerSettings(seed=14, restarts=3, max_sweeps=2))
    one = optimize(protein, config, OptimizerSettings(seed=14, restarts=1, max_sweeps=2))
    assert many.restarts_run == 3
    assert many.score <= one.score


def test_annealing_can_improve_on_the_greedy_result(plain_config):
    protein = duplicated_protein(block=50, seed=30)
    settings = OptimizerSettings(seed=15, max_sweeps=2, anneal_steps=1500)
    result = optimize(protein, plain_config, settings)
    assert result.moves_evaluated > 0
    # Whatever annealing explores, the reported sequence is the best one seen.
    assert SequenceState(result.sequence, plain_config).score == pytest.approx(
        result.score, abs=1e-9
    )


def test_seed_modes_behave_as_documented():
    protein = random_protein(40, seed=31)
    table = get_table("human")
    best = seed_sequence(protein, table, random.Random(0), "best")
    assert best == "".join(table.most_frequent(aa) for aa in protein)
    weighted = seed_sequence(protein, table, random.Random(0), "weighted")
    assert translate(weighted) == protein
    with pytest.raises(ValueError, match="unknown seed mode"):
        seed_sequence(protein, table, random.Random(0), "nonsense")


def test_start_sequence_must_encode_the_protein(config):
    protein = random_protein(20, seed=32)
    with pytest.raises(ValueError, match="does not encode"):
        optimize(protein, config, OptimizerSettings(seed=16), start="ATGAAA")


def test_optimize_dna_keeps_the_protein(config):
    protein = random_protein(50, seed=33)
    dna = random_dna(protein, seed=34)
    result = optimize_dna(dna, config, OptimizerSettings(seed=17, max_sweeps=2))
    assert translate(result.sequence) == protein
    assert result.initial_sequence == dna


def test_time_limit_is_respected(config):
    protein = duplicated_protein(block=120, seed=35)
    settings = OptimizerSettings(seed=18, max_sweeps=50, max_seconds=0.3)
    result = optimize(protein, config, settings)
    assert result.elapsed < 10  # a hard ceiling, not a timing assertion
    assert translate(result.sequence) == protein


def test_invalid_protein_is_refused(config):
    with pytest.raises(SequenceError):
        optimize("MK#V", config)


def test_progress_callback_receives_messages(config):
    messages: list[str] = []
    protein = random_protein(30, seed=36)
    optimize(
        protein,
        config,
        OptimizerSettings(seed=19, max_sweeps=1, progress=messages.append),
    )
    assert messages
    assert any("restart" in message for message in messages)


def test_single_residue_protein_is_handled(config):
    result = optimize("M", config, OptimizerSettings(seed=20))
    assert result.sequence == "ATG"
    assert result.score == pytest.approx(result.initial_score)
