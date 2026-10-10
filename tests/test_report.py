"""Tests for the before/after reporting."""

from __future__ import annotations

import json

from codonopt import OptimizerSettings, format_report, measure, optimize, report_dict

from conftest import duplicated_pair


def test_measure_covers_every_constraint(config):
    _, dna = duplicated_pair(block=40, seed=40)
    summary = measure(dna, config)
    for key in (
        "length",
        "gc_fraction",
        "cai",
        "redundant_kmers",
        "hairpin_pairs",
        "longest_homopolymer",
        "forbidden_motif_hits",
        "splice_backend",
        "splice_donor_sites",
    ):
        assert key in summary
    assert summary["length"] == len(dna)
    assert summary["redundant_kmers"] > 0


def test_measure_omits_splice_keys_when_scoring_is_off(plain_config):
    summary = measure("ATGAAATGGGTT", plain_config)
    assert "splice_backend" not in summary
    assert "splice_donor_sites" not in summary


def test_longest_homopolymer_reports_runs_under_the_limit(plain_config):
    # Three A's is under the default limit of six but should still be shown.
    summary = measure("ATGAAACCC", plain_config.with_(max_homopolymer=6))
    assert summary["longest_homopolymer"] == 3
    assert summary["homopolymer_runs_over_limit"] == 0
    assert summary["homopolymer_excess_bases"] == 0


def test_report_dict_is_json_serialisable(config):
    protein, dna = duplicated_pair(block=30, seed=41)
    result = optimize(
        protein, config, OptimizerSettings(seed=42, max_sweeps=2), start=dna
    )
    payload = report_dict(result, config)
    text = json.dumps(payload)
    assert json.loads(text)["sequence"] == result.sequence
    assert payload["initial"]["redundant_kmers"] >= payload["final"]["redundant_kmers"]
    assert payload["search"]["restarts_run"] == 1
    assert set(payload["score"]["final_terms"]) == set(result.breakdown.terms)


def test_format_report_shows_before_and_after(config):
    protein, dna = duplicated_pair(block=30, seed=43)
    result = optimize(
        protein, config, OptimizerSettings(seed=44, max_sweeps=2), start=dna
    )
    text = format_report(result, config)
    assert "before" in text and "after" in text
    assert "Redundant k-mers" in text
    assert "Splice donor sites" in text
    assert f"{len(protein)} residues" in text
