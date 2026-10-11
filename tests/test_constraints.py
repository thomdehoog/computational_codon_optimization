"""Regression cases from the algorithm review, including independent oracles."""

import itertools
import json
import math
import random
import importlib

import pytest

from codonopt import (DesignConfig, OptimizerSettings, Weights, Motif, SequenceState,
                      constraint_violations, optimize, optimize_dna, get_table,
                      resolve_motifs, report_dict, measure, translate)
from codonopt.cli import main
from codonopt.sequence import reverse_complement
from codonopt.splice import MotifScorer, get_scorer
from codonopt.tables import CodonUsage


def only(**terms):
    return Weights(**{name: terms.get(name, 0) for name in vars(Weights())})


def test_coordinated_repair_reaches_exhaustive_feasible_optimum():
    config = DesignConfig(table=get_table("uniform"), weights=only(motif=1),
                          motifs=tuple(resolve_motifs(["GAATTC", "GAGTTC", "GAATTT"])))
    result = optimize_dna("GAATTC", config, OptimizerSettings(seed=1, max_seconds=None))
    candidates = [a+b for a, b in itertools.product(("GAA", "GAG"), ("TTT", "TTC"))]
    feasible = [s for s in candidates if not constraint_violations(s, config)]
    assert result.valid
    assert result.sequence == feasible[0] == "GAGTTT"
    assert result.score == min(SequenceState(s, config).score for s in candidates)
    assert result.repair_evaluated > 0


def test_validity_outranks_score_and_zero_weight_does_not_disable_a_requirement():
    config = DesignConfig(motifs=tuple(resolve_motifs(["TTC"])), weights=only(usage=1000))
    result = optimize_dna("TTC", config, OptimizerSettings(seed=1))
    assert result.valid and result.sequence == "TTT"
    assert result.score > result.initial_score
    assert not result.unresolved_constraints


def test_unavoidable_constraint_is_unresolved_not_declared_impossible():
    config = DesignConfig(motifs=(Motif("start", "ATG"),))
    result = optimize("MW", config, OptimizerSettings(seed=1))
    assert not result.valid
    assert result.unresolved_constraints[0]["criterion"] == "motif"
    assert result.unresolved_constraints[0]["count"] == 1
    assert result.stop_reason == "completed"
    assert "not proof" in result.notes[0]


def test_cli_marks_unresolved_output_even_when_quiet(capsys):
    assert main(["optimize", "--protein", "MW", "--avoid", "ATG", "--quiet"]) == 2
    output = capsys.readouterr()
    assert "valid=false" in output.out
    assert "unresolved" in output.err


def test_gc_has_the_same_per_base_tradeoff_across_lengths():
    config = DesignConfig(weights=only(gc_global=50, usage=.1))
    for length in (100, 1000):
        result = optimize("C" * length, config, OptimizerSettings(seed=1, max_seconds=None))
        assert measure(result.sequence, config)["gc_fraction"] == .5


def test_flanks_are_scored_and_immutable_with_an_unaligned_cds():
    config = DesignConfig(left_flank="GGT", right_flank="A",
                          motifs=tuple(resolve_motifs(["BsaI"])))
    result = optimize_dna("CTC", config, OptimizerSettings(seed=1))
    assert result.valid
    assert translate(result.sequence) == "L"
    assert result.initial_sequence == "CTC"
    assert measure(result.initial_sequence, config)["forbidden_motif_hits"] == 1
    assert measure(result.sequence, config)["forbidden_motif_hits"] == 0
    assert report_dict(result, config)["context_sequence"] == "GGT" + result.sequence + "A"


def test_fixed_only_violations_are_reported_without_editing_flanks():
    config = DesignConfig(left_flank="GAATTC", motifs=tuple(resolve_motifs(["EcoRI"])))
    result = optimize("M", config, OptimizerSettings(seed=1))
    assert result.sequence == "ATG" and not result.valid
    assert result.unresolved_constraints[0]["regions"][0]["start"] == 0


def test_protected_codons_survive_restarts_and_annealing():
    config = DesignConfig(protected_codons=(0, 2))
    result = optimize_dna("TTTCTATTT", config,
                          OptimizerSettings(seed=0, restarts=3, anneal_steps=100))
    assert result.sequence[:3] == result.sequence[-3:] == "TTT"
    assert result.initial_sequence == "TTTCTATTT"
    with pytest.raises(ValueError, match="require a starting"):
        optimize("FLF", config)


def test_original_baseline_is_preserved_when_a_later_restart_wins():
    protein = "M" + "".join(random.Random(41).choices("ACDEFGHIKLMNPQRSTVWY", k=59))
    dna = "".join(get_table("human").most_frequent(a) for a in protein)
    config = DesignConfig()
    result = optimize_dna(dna, config, OptimizerSettings(seed=0, restarts=3, max_sweeps=3))
    assert result.initial_sequence == dna
    assert result.initial_score == pytest.approx(SequenceState(dna, config).score)
    assert report_dict(result, config)["initial_sequence"] == dna
    assert translate(result.winning_seed_sequence) == protein


def test_longest_hairpin_comes_from_the_complete_scan():
    stem = "".join(random.Random(0).choices("ACGT", k=28))
    unit = "ACGTTGCACC" + "TTTTT" + reverse_complement("ACGTTGCACC") + "AAA"
    sequence = unit*12 + "AAA" + stem + "AAA" + reverse_complement(stem)
    sequence += "A" * (-len(sequence) % 3)
    assert measure(sequence, DesignConfig())["longest_hairpin_stem"] >= 28


def test_annealing_schedule_cools_on_noop_draws(monkeypatch):
    module = importlib.import_module("codonopt.optimize")
    state = SequenceState("AAA", DesignConfig(weights=only(motif=1), required=(),
                          motifs=(Motif("uphill", "AAG"),)))
    class Rng:
        draws = 0
        def choice(self, options):
            if options == [0]:
                return 0
            self.draws += 1
            return "AAA" if self.draws % 2 else "AAG"
        def random(self):
            return 1
    exponents = []
    original_exp = math.exp
    def capture(x):
        if x < -1:
            exponents.append(x)
        return original_exp(x)
    monkeypatch.setattr(module.math, "exp", capture)
    module._anneal(state, Rng(), OptimizerSettings(anneal_steps=10, start_temperature=2,
                   end_temperature=.05), lambda: False)
    assert exponents[-1] == pytest.approx(-20)


def test_feasible_annealing_never_reintroduces_required_motifs():
    config = DesignConfig(motifs=(Motif("bad", "AAA"),), weights=only())
    result = optimize_dna("AAGAAG", config,
                          OptimizerSettings(seed=8, anneal_steps=500, start_temperature=1e4))
    assert result.valid and not constraint_violations(result.sequence, config)


def test_splice_checks_expressed_strand_by_default_and_auto_is_host_aware():
    config = DesignConfig(splice_scorer=MotifScorer(), weights=only(splice=1))
    assert SequenceState("ACTTACCTG", config).score == 0
    assert SequenceState("ACTTACCTG", config.with_(splice_strand="both")).score == 3
    assert get_scorer("auto", host="ecoli") is None
    assert get_scorer("auto", host="custom") is None
    assert get_scorer("motif", host="ecoli") is not None  # explicit override


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_core_values_are_rejected(value):
    for make in (lambda: Weights(usage=value),
                 lambda: OptimizerSettings(start_temperature=value),
                 lambda: OptimizerSettings(max_seconds=value),
                 lambda: DesignConfig(donor_cutoff=value),
                 lambda: CodonUsage.from_frequencies("bad", {"TTT": value})):
        with pytest.raises(ValueError):
            make()


def test_motifs_are_canonicalized_and_conflicting_names_are_rejected():
    config = DesignConfig(motifs=(Motif("one", "GAATTC"), Motif("two", "GAATTC")))
    state = SequenceState("GAATTC", config)
    assert state.breakdown().total == pytest.approx(state.score)
    assert len(config.motifs) == 1
    with pytest.raises(ValueError, match="different patterns"):
        resolve_motifs(["same=GAATTC", "same=GGATCC"])
    with pytest.raises(ValueError, match="different patterns"):
        DesignConfig(motifs=(Motif("same", "GAATTC"), Motif("same", "GGATCC")))


def test_excluded_codons_are_hard_requirements_even_with_one_alternative():
    config = DesignConfig(excluded_codons=("TTT",), weights=only())
    result = optimize_dna("TTT", config, OptimizerSettings(seed=1, restarts=2))
    assert result.sequence == "TTC" and result.valid
    with pytest.raises(ValueError, match="all codons"):
        optimize("F", config.with_(excluded_codons=("TTT", "TTC")))
    locked = optimize_dna("TTT", config.with_(protected_codons=(0,)),
                          OptimizerSettings(seed=1, restarts=2))
    assert not locked.valid and locked.sequence == "TTT"


def test_saved_record_can_replay_a_generated_seed():
    config = DesignConfig(motifs=tuple(resolve_motifs(["EcoRI"])), left_flank="G")
    settings = OptimizerSettings(max_seconds=None, restarts=2, anneal_steps=100)
    result = optimize("MELF", config, settings)
    payload = report_dict(result, config)
    json.dumps(payload, allow_nan=False)
    # Reconstruct from JSON bytes, not the original in-memory configuration.
    saved = json.loads(json.dumps(payload, allow_nan=False))
    values = saved["configuration"]
    table = values.pop("table")
    table.pop("sha256")
    backend = values.pop("splice_scorer")
    values["weights"] = Weights(**values["weights"])
    values["motifs"] = tuple(Motif(m["name"], m["pattern"]) for m in values["motifs"])
    restored = DesignConfig(**values, table=CodonUsage(**table),
                            splice_scorer=get_scorer(backend["name"]))
    replay = optimize("MELF", restored, OptimizerSettings(**saved["search"]["settings"]))
    assert replay.sequence == result.sequence
    assert payload["configuration"]["table"]["frequencies"] == config.table.frequencies
    assert len(payload["configuration"]["table"]["sha256"]) == 64
    assert payload["configuration"]["splice_scorer"]["name"] == "none"


def test_incremental_constraints_with_flanks_match_independent_validation():
    config = DesignConfig(left_flank="AA", right_flank="CGTGC", repeat_k=3,
                          hairpin_stem=2, hairpin_min_loop=0, hairpin_max_loop=10,
                          max_homopolymer=2, gc_window=5, gc_step=1,
                          splice_scorer=MotifScorer(), splice_strand="both",
                          motifs=(Motif("pattern", "AA"),),
                          required=("motif", "gc_window", "repeat", "hairpin", "homopolymer", "splice"),
                          weights=only())
    state = SequenceState("AAAACGGCTTTC", config)
    rng = random.Random(57)
    for _ in range(150):
        i = rng.randrange(state.codon_count)
        state.set_codon(i, rng.choice(state.alternatives(i)))
        counts = sum(v["count"] for v in constraint_violations(state.sequence, config))
        assert state.violations == counts
        assert state.score == pytest.approx(state.audit())
        assert state.breakdown().total == pytest.approx(state.score)
        assert state.context_sequence[:2] == "AA"
        assert state.context_sequence[-5:] == "CGTGC"


def test_budget_and_cancellation_keep_a_consistent_result():
    config = DesignConfig(motifs=tuple(resolve_motifs(["GAATTC", "GAGTTC", "GAATTT"])),
                          weights=only(motif=1))
    limited = optimize_dna("GAATTC", config, OptimizerSettings(seed=1, repair_candidates=1))
    assert limited.repair_evaluated <= 1
    assert limited.valid == (not constraint_violations(limited.sequence, config))
    calls = 0
    def pulse():
        nonlocal calls
        calls += 1
        return calls == 1
    cancelled = optimize_dna("GAATTC", config, OptimizerSettings(seed=1, should_stop=pulse))
    assert cancelled.cancelled and cancelled.stop_reason == "cancelled"
    no_search = optimize("M", settings=OptimizerSettings(max_sweeps=0, should_stop=lambda: True))
    assert no_search.cancelled and no_search.valid


def test_cli_add_stop_keeps_original_dna_as_baseline(tmp_path):
    path = tmp_path / "result.json"
    assert main(["optimize", "--dna", "TTTCTATTT", "--add-stop", "--sweeps", "0",
                 "--json", str(path), "--quiet", "--seed", "1"]) == 0
    payload = json.loads(path.read_text())
    assert payload["initial_sequence"].startswith("TTTCTATTT")


@pytest.mark.parametrize("criterion", ["gc_window", "repeat", "hairpin", "homopolymer", "splice"])
def test_required_criteria_work_when_the_soft_weight_is_zero(criterion):
    config = DesignConfig(required=(criterion,), weights=only(), gc_window=6,
                          gc_range=(.4, .6), repeat_k=3, hairpin_stem=2,
                          hairpin_min_loop=0, hairpin_max_loop=10,
                          max_homopolymer=2, splice_scorer=MotifScorer())
    state = SequenceState("CAGGTAAGTAAAAAATTT", config)
    assert state.violations == sum(v["count"] for v in constraint_violations(state.sequence, config))
