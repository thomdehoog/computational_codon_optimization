"""Tests for the command line interface."""

from __future__ import annotations

import json

import pytest

from codonopt.cli import main
from codonopt.fasta import read_fasta
from codonopt.sequence import translate

PROTEIN = (
    "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEENFKALVLIAFAQYLQQCPF"
    "DEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEP"
)
FAST = ["--splice-backend", "none", "--seed", "1", "--sweeps", "2"]


def test_optimize_writes_fasta_and_json(tmp_path, capsys):
    protein_path = tmp_path / "p.faa"
    protein_path.write_text(f">albumin fragment\n{PROTEIN}\n")
    out = tmp_path / "out.fna"
    report = tmp_path / "out.json"

    code = main(
        [
            "optimize",
            "--protein-fasta", str(protein_path),
            "--out", str(out),
            "--json", str(report),
            "--avoid", "golden_gate,EcoRI",
            *FAST,
        ]
    )
    assert code == 0

    record = read_fasta(out)[0]
    assert record.name == "albumin"
    assert translate(record.sequence) == PROTEIN

    payload = json.loads(report.read_text())
    assert payload["final"]["forbidden_motif_hits"] == 0
    assert payload["score"]["final"] <= payload["score"]["initial"]
    assert payload["sequence"] == record.sequence
    assert payload["codon_usage_table"] == "human"


def test_optimize_writes_to_stdout_when_no_out_given(capsys):
    code = main(["optimize", "--protein", "MKWV", "--quiet", *FAST])
    assert code == 0
    captured = capsys.readouterr()
    assert captured.out.startswith(">")
    assert translate(captured.out.splitlines()[1]) == "MKWV"


def test_quiet_suppresses_the_report(capsys):
    main(["optimize", "--protein", "MKWV", "--quiet", *FAST])
    assert capsys.readouterr().err == ""
    main(["optimize", "--protein", "MKWV", *FAST])
    assert "Objective" in capsys.readouterr().err


def test_optimize_from_dna_input(tmp_path, capsys):
    code = main(["optimize", "--dna", "ATGAAATGGGTT", "--quiet", *FAST])
    assert code == 0
    sequence = capsys.readouterr().out.splitlines()[1]
    assert translate(sequence) == "MKWV"


def test_add_stop_appends_a_stop_codon(capsys):
    main(["optimize", "--protein", "MKWV", "--add-stop", "--quiet", *FAST])
    sequence = capsys.readouterr().out.splitlines()[1]
    assert translate(sequence) == "MKWV*"


def test_name_option_overrides_the_record_name(capsys):
    main(["optimize", "--protein", "MKWV", "--name", "construct", "--quiet", *FAST])
    assert capsys.readouterr().out.startswith(">construct")


def test_weight_overrides_are_applied(capsys):
    # With only codon usage weighted, the result is the most frequent codon
    # at every position.
    main(
        [
            "optimize", "--protein", "LLLL", "--quiet", "--table", "human",
            "--splice-backend", "none", "--seed", "1",
            "--weight", "gc_global=0,gc_window=0,repeat=0,hairpin=0",
            "--weight", "homopolymer=0", "--weight", "usage=1",
        ]
    )
    assert capsys.readouterr().out.splitlines()[1] == "CTG" * 4


def test_report_subcommand_prints_measurements(capsys):
    code = main(["report", "--dna", "ATGAAATGGGTTACC", "--splice-backend", "none"])
    assert code == 0
    output = capsys.readouterr().out
    assert "gc_fraction" in output
    assert "protein_length" in output


def test_report_subcommand_writes_json(tmp_path):
    report = tmp_path / "m.json"
    main(
        [
            "report", "--dna", "ATGAAATGGGTTACC",
            "--json", str(report), "--splice-backend", "none",
        ]
    )
    payload = json.loads(report.read_text())
    assert payload["length"] == 15
    assert payload["protein_length"] == 5


def test_tables_and_motifs_subcommands(capsys):
    assert main(["tables"]) == 0
    assert "human" in capsys.readouterr().out
    assert main(["motifs"]) == 0
    output = capsys.readouterr().out
    assert "golden_gate" in output
    assert "EcoRI" in output


def test_table_file_option(tmp_path, capsys):
    table = tmp_path / "usage.json"
    table.write_text(json.dumps({"CTT": 1.0}))
    main(
        [
            "optimize", "--protein", "LL", "--quiet", "--table-file", str(table),
            "--splice-backend", "none", "--seed", "1",
            "--weight", "gc_global=0,gc_window=0,usage=1",
        ]
    )
    assert capsys.readouterr().out.splitlines()[1] == "CTTCTT"


@pytest.mark.parametrize(
    "argv, message",
    [
        (["optimize", "--protein", "MKZ"], "unsupported residues"),
        (["optimize", "--dna", "ATGA"], "multiple of three"),
        (["optimize", "--protein", "MK", "--gc-range", "30", "70"], "fractions"),
        (["optimize", "--protein", "MK", "--avoid", "NotAnEnzyme"], "cannot resolve"),
        (["optimize", "--protein", "MK", "--weight", "nonsense=1"], "unknown weight"),
        (["optimize", "--protein", "MK", "--weight", "usage"], "NAME=VALUE"),
        (["optimize", "--protein", "MK", "--weight", "usage=high"], "not a number"),
        (["report", "--dna-fasta", "/nonexistent/path.fna"], ""),
        (["optimize", "--protein", ""], "empty protein sequence"),
        (["optimize", "--dna", ""], "empty nucleotide sequence"),
    ],
)
def test_errors_are_reported_cleanly(argv, message, capsys):
    assert main([*argv, *(["--splice-backend", "none"] if argv[0] == "optimize" else [])]) == 1
    captured = capsys.readouterr()
    assert captured.err.startswith("codonopt: error:")
    assert message in captured.err


def test_multi_record_fasta_is_refused(tmp_path, capsys):
    path = tmp_path / "two.faa"
    path.write_text(">a\nMKV\n>b\nMKW\n")
    assert main(["optimize", "--protein-fasta", str(path), "--splice-backend", "none"]) == 1
    assert "one at a time" in capsys.readouterr().err


def test_input_flags_are_mutually_exclusive():
    with pytest.raises(SystemExit):
        main(["optimize", "--protein", "MK", "--dna", "ATG"])


def test_missing_input_is_a_usage_error():
    with pytest.raises(SystemExit):
        main(["optimize"])


def test_report_does_not_offer_weight_overrides(capsys):
    """`report` measures a sequence, so penalty weights do not apply to it."""
    with pytest.raises(SystemExit):
        main(["report", "--dna", "ATGAAA", "--weight", "usage=1"])


def test_empty_sequence_is_a_clean_error_not_a_traceback(capsys):
    assert main(["report", "--dna", "", "--splice-backend", "none"]) == 1
    assert "empty nucleotide sequence" in capsys.readouterr().err
