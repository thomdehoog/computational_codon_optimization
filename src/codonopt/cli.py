"""Command-line interface.

    codonopt optimize --protein-fasta gene.faa --table human --out gene.fna
    codonopt report --dna-fasta gene.fna
    codonopt gui
    codonopt tables
    codonopt motifs
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import fields, replace
from pathlib import Path
from typing import Sequence

from . import __version__
from .config import DesignConfig, OptimizerSettings, Weights
from .constraints import unresolved_is_all_fixed_context
from .fasta import Record, format_fasta, read_fasta
from .motifs import MOTIF_GROUPS, REGULATORY_MOTIFS, RESTRICTION_SITES, resolve_motifs
from .optimize import optimize
from .report import format_report, measure, report_dict
from .sequence import SequenceError, translate, validate_dna, validate_protein
from .splice import get_scorer, maxent_available
from .tables import available_tables, get_table, load_table

__all__ = ["main", "build_parser"]


def build_parser() -> argparse.ArgumentParser:
    """Return the fully configured argument parser."""
    parser = argparse.ArgumentParser(
        prog="codonopt",
        description="Codon-optimize a coding sequence for synthesis and expression.",
    )
    parser.add_argument("--version", action="version", version=f"codonopt {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    optimize_parser = subparsers.add_parser(
        "optimize", help="design a coding sequence for a protein"
    )
    _add_input_arguments(optimize_parser)
    _add_constraint_arguments(optimize_parser)
    _add_search_arguments(optimize_parser)
    optimize_parser.add_argument(
        "--out",
        metavar="PATH",
        help="write the optimized sequence here as FASTA (default: stdout)",
    )
    optimize_parser.add_argument(
        "--json", metavar="PATH", help="write a machine-readable report here"
    )
    optimize_parser.add_argument(
        "--name",
        help="FASTA record name for the output; the default is the input "
        "record's name, or 'protein'/'sequence' for a sequence given inline",
    )
    optimize_parser.add_argument(
        "--quiet", action="store_true", help="suppress the progress and summary report"
    )
    optimize_parser.set_defaults(handler=_run_optimize)

    report_parser = subparsers.add_parser(
        "report", help="measure an existing coding sequence"
    )
    _add_input_arguments(report_parser, protein_allowed=False)
    _add_constraint_arguments(report_parser, weights=False)
    report_parser.add_argument(
        "--json", metavar="PATH", help="write the measurements here as JSON"
    )
    report_parser.set_defaults(handler=_run_report)

    gui_parser = subparsers.add_parser(
        "gui",
        help="open the graphical interface (a native window with pywebview, "
        "otherwise your browser)",
    )
    gui_parser.add_argument(
        "--port", type=int, default=0, metavar="N", help="port to listen on (default: any free one)"
    )
    gui_parser.add_argument(
        "--no-window",
        action="store_true",
        help="serve the interface without opening a native window",
    )
    gui_parser.add_argument(
        "--no-browser",
        action="store_true",
        help="do not open the default browser (use with --no-window)",
    )
    gui_parser.add_argument(
        "--verbose", action="store_true", help="log each request to stderr"
    )
    gui_parser.set_defaults(handler=_run_gui)

    tables_parser = subparsers.add_parser("tables", help="list codon usage tables")
    tables_parser.set_defaults(handler=_run_tables)

    motifs_parser = subparsers.add_parser("motifs", help="list known motifs and groups")
    motifs_parser.set_defaults(handler=_run_motifs)
    return parser


def _add_input_arguments(
    parser: argparse.ArgumentParser, *, protein_allowed: bool = True
) -> None:
    group = parser.add_mutually_exclusive_group(required=True)
    if protein_allowed:
        group.add_argument("--protein", metavar="SEQ", help="protein sequence")
        group.add_argument(
            "--protein-fasta", metavar="PATH", help="FASTA file holding the protein"
        )
    group.add_argument(
        "--dna",
        metavar="SEQ",
        help="in-frame coding sequence"
        + (" to improve" if protein_allowed else " to measure"),
    )
    group.add_argument(
        "--dna-fasta", metavar="PATH", help="FASTA file holding the coding sequence"
    )
    if protein_allowed:
        parser.add_argument(
            "--add-stop",
            action="store_true",
            help="append a stop codon if the protein does not end in '*'",
        )


def _add_constraint_arguments(
    parser: argparse.ArgumentParser, *, weights: bool = True
) -> None:
    group = parser.add_argument_group("design constraints")
    group.add_argument(
        "--table",
        default="human",
        help="bundled codon usage table: " + ", ".join(available_tables()),
    )
    group.add_argument(
        "--table-file",
        metavar="PATH",
        help="codon usage table as JSON, CSV or TSV, overriding --table",
    )
    group.add_argument(
        "--gc-target", type=float, default=0.50, metavar="F", help="target GC fraction"
    )
    group.add_argument(
        "--gc-range",
        type=float,
        nargs=2,
        default=(0.30, 0.70),
        metavar=("LOW", "HIGH"),
        help="acceptable GC fraction per window",
    )
    group.add_argument("--gc-window", type=int, default=30, metavar="N")
    group.add_argument("--gc-step", type=int, default=3, metavar="N")
    group.add_argument(
        "--repeat-k",
        type=int,
        default=12,
        metavar="N",
        help="k-mer size treated as a direct repeat",
    )
    group.add_argument(
        "--hairpin-stem",
        type=int,
        default=10,
        metavar="N",
        help="stem length treated as a hairpin arm",
    )
    group.add_argument(
        "--hairpin-loop",
        type=int,
        nargs=2,
        default=(3, 60),
        metavar=("MIN", "MAX"),
        help="loop lengths counted as a hairpin",
    )
    group.add_argument(
        "--max-homopolymer",
        type=int,
        default=6,
        metavar="N",
        help="longest acceptable single-base run",
    )
    group.add_argument(
        "--avoid",
        action="append",
        default=[],
        metavar="SPEC",
        help="motifs to exclude: a group, an enzyme name, an IUPAC pattern, or "
        "NAME=PATTERN; repeat or comma-separate (see 'codonopt motifs')",
    )
    group.add_argument(
        "--splice-backend",
        default="none",
        choices=("auto", "maxent", "motif", "none"),
        help="candidate splice signal screening (default: none; auto is human-only)",
    )
    group.add_argument(
        "--donor-cutoff",
        type=float,
        metavar="F",
        help="highest acceptable donor score (default: the backend's)",
    )
    group.add_argument(
        "--acceptor-cutoff",
        type=float,
        metavar="F",
        help="highest acceptable acceptor score (default: the backend's)",
    )
    group.add_argument("--splice-strand", choices=("forward", "reverse", "both"), default="forward")
    group.add_argument("--require", action="append", default=[],
                       choices=("motif", "gc_window", "repeat", "hairpin", "homopolymer", "splice"),
                       help="make a criterion mandatory; repeat to require several")
    group.add_argument("--soft-motifs", action="store_true", help="treat avoided motifs as preferences instead of requirements")
    group.add_argument("--left-flank", default="", help="immutable upstream DNA used when scoring the construct")
    group.add_argument("--right-flank", default="", help="immutable downstream DNA used when scoring the construct")
    group.add_argument("--protect-codon", type=int, action="append", default=[], metavar="INDEX",
                       help="protect a zero-based codon index in starting DNA; repeat as needed")
    group.add_argument("--exclude-codon", action="append", default=[], metavar="CODON",
                       help="forbid a codon independently of its frequency")
    if weights:
        group.add_argument(
            "--weight",
            action="append",
            default=[],
            metavar="NAME=VALUE",
            help="override a penalty weight: "
            + ", ".join(f.name for f in fields(Weights)),
        )


def _add_search_arguments(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("search effort")
    group.add_argument("--restarts", type=int, default=1, metavar="N")
    group.add_argument(
        "--sweeps", type=int, default=6, metavar="N", help="greedy passes per restart"
    )
    group.add_argument(
        "--anneal",
        type=int,
        default=0,
        metavar="N",
        help="simulated annealing steps after the greedy passes",
    )
    group.add_argument("--temp-start", type=float, default=2.0, metavar="F")
    group.add_argument("--temp-end", type=float, default=0.05, metavar="F")
    group.add_argument(
        "--seconds",
        type=float,
        default=120.0,
        metavar="F",
        help="wall-clock budget; 0 means no limit",
    )
    group.add_argument(
        "--seed", type=int, metavar="N", help="random seed, for reproducible runs"
    )
    group.add_argument(
        "--seed-mode", default="weighted", choices=("weighted", "best", "random")
    )
    group.add_argument("--repair-candidates", type=int, default=2000)
    group.add_argument("--repair-max-codons", type=int, default=3, choices=(2, 3, 4))


class _UsageError(Exception):
    """Raised for input the user can fix."""


def _build_config(args: argparse.Namespace) -> DesignConfig:
    table = (
        load_table(args.table_file) if args.table_file else get_table(args.table)
    )
    try:
        scorer = get_scorer(args.splice_backend, host=table.name)
    except RuntimeError as exc:
        raise _UsageError(str(exc)) from None

    specs: list[str] = []
    for entry in args.avoid:
        specs.extend(part for part in entry.split(",") if part.strip())
    motifs = tuple(resolve_motifs(specs))

    config = DesignConfig(
        table=table,
        gc_target=args.gc_target,
        gc_range=(args.gc_range[0], args.gc_range[1]),
        gc_window=args.gc_window,
        gc_step=args.gc_step,
        repeat_k=args.repeat_k,
        hairpin_stem=args.hairpin_stem,
        hairpin_min_loop=args.hairpin_loop[0],
        hairpin_max_loop=args.hairpin_loop[1],
        max_homopolymer=args.max_homopolymer,
        motifs=motifs,
        splice_scorer=scorer,
        donor_cutoff=args.donor_cutoff,
        acceptor_cutoff=args.acceptor_cutoff,
        weights=_build_weights(getattr(args, "weight", [])),
        required=tuple(args.require) + (() if args.soft_motifs else ("motif",)),
        splice_strand=args.splice_strand,
        left_flank=args.left_flank,
        right_flank=args.right_flank,
        protected_codons=tuple(args.protect_codon),
        excluded_codons=tuple(args.exclude_codon),
    )
    return config


def _build_weights(overrides: list[str]) -> Weights:
    known = {field.name for field in fields(Weights)}
    changes: dict[str, float] = {}
    for item in overrides:
        for part in item.split(","):
            part = part.strip()
            if not part:
                continue
            name, separator, raw = part.partition("=")
            if not separator:
                raise _UsageError(f"--weight expects NAME=VALUE, got {part!r}")
            name = name.strip()
            if name not in known:
                raise _UsageError(
                    f"unknown weight {name!r}; choose from {', '.join(sorted(known))}"
                )
            try:
                changes[name] = float(raw)
            except ValueError:
                raise _UsageError(f"--weight {name}: {raw!r} is not a number") from None
    return replace(Weights(), **changes)


def _read_input(args: argparse.Namespace) -> tuple[str, str, bool]:
    """Return ``(sequence, record name, is_dna)`` from the chosen input flag.

    Dispatch is on which flag was given rather than on whether its value is
    non-empty, so ``--protein ""`` reaches sequence validation and reports an
    empty sequence instead of falling through to another branch.
    """
    if getattr(args, "protein", None) is not None:
        return args.protein, "protein", False
    if getattr(args, "protein_fasta", None) is not None:
        record = _single_record(args.protein_fasta)
        return record.sequence, record.name, False
    if args.dna is not None:
        return args.dna, "sequence", True
    record = _single_record(args.dna_fasta)
    return record.sequence, record.name, True


def _single_record(path: str) -> Record:
    records = read_fasta(path)
    if len(records) > 1:
        raise _UsageError(
            f"{path}: holds {len(records)} records; this tool handles one at a time"
        )
    return records[0]


def _run_optimize(args: argparse.Namespace) -> int:
    config = _build_config(args)
    raw, name, is_dna = _read_input(args)

    start = None
    if is_dna:
        start = validate_dna(raw, in_frame=True)
        protein = translate(start)
    else:
        protein = validate_protein(raw)
    if args.add_stop and not protein.endswith("*"):
        protein += "*"
        if start is not None:
            options = [c for c in config.table.codons("*") if c not in config.excluded_codons]
            if not options:
                raise _UsageError("all stop codons are excluded")
            start += max(options, key=config.table.frequency)

    log = (lambda message: None) if args.quiet else _stderr
    settings = OptimizerSettings(
        restarts=args.restarts,
        max_sweeps=args.sweeps,
        anneal_steps=args.anneal,
        start_temperature=args.temp_start,
        end_temperature=args.temp_end,
        max_seconds=None if args.seconds == 0 else args.seconds,
        seed=args.seed,
        seed_mode=args.seed_mode,
        progress=log,
        repair_candidates=args.repair_candidates,
        repair_max_codons=args.repair_max_codons,
    )
    if config.splice_scorer is not None and not args.quiet:
        log(f"splice backend: {config.splice_scorer.name}")

    result = optimize(protein, config, settings, start=start)

    record = Record(
        name=args.name or name,
        description=(args.name or name)
            + f" codonopt table={config.table.name} score={result.score:.2f} valid={str(result.valid).lower()}",
        sequence=result.sequence,
    )
    text = format_fasta(record)
    if args.out:
        Path(args.out).write_text(text)
        log(f"wrote {args.out}")
    else:
        sys.stdout.write(text)
    if args.json:
        Path(args.json).write_text(
            json.dumps(report_dict(result, config), indent=2, allow_nan=False) + "\n"
        )
        log(f"wrote {args.json}")
    if not args.quiet:
        log("")
        log(format_report(result, config))
    if not result.valid:
        if unresolved_is_all_fixed_context(result.unresolved_constraints):
            _stderr(
                "codonopt: the coding sequence meets every requirement, but the fixed "
                "flanks you supplied do not; no codon choice can change that"
            )
        else:
            _stderr("codonopt: unresolved mandatory constraints; output is a candidate, not a valid design")
        return 2
    return 0


def _run_report(args: argparse.Namespace) -> int:
    config = _build_config(args)
    raw, name, _ = _read_input(args)
    sequence = validate_dna(raw, in_frame=True)
    summary = measure(sequence, config)
    summary["protein_length"] = len(translate(sequence))

    if args.json:
        Path(args.json).write_text(json.dumps(summary, indent=2) + "\n")
    width = max(len(key) for key in summary)
    print(f"{name}:")
    for key, value in summary.items():
        if isinstance(value, list):
            value = ", ".join(str(item) for item in value) or "-"
        print(f"  {key:<{width}}  {value}")
    return 0


def _run_gui(args: argparse.Namespace) -> int:
    from .server import run_gui

    return run_gui(
        args.port,
        window=not args.no_window,
        open_browser=not args.no_browser,
        verbose=args.verbose,
    )


def _run_tables(args: argparse.Namespace) -> int:
    for name in available_tables():
        table = get_table(name)
        print(f"{name:<10}{table.description}")
    return 0


def _run_motifs(args: argparse.Namespace) -> int:
    print("groups:")
    for name, members in MOTIF_GROUPS.items():
        print(f"  {name:<16}{', '.join(members) if members else '(empty)'}")
    print("\nrestriction sites:")
    for name, pattern in sorted(RESTRICTION_SITES.items()):
        print(f"  {name:<10}{pattern}")
    print("\nregulatory motifs:")
    for name, pattern in sorted(REGULATORY_MOTIFS.items()):
        print(f"  {name:<18}{pattern}")
    print(
        "\nAny IUPAC pattern or NAME=PATTERN also works, e.g. "
        "--avoid golden_gate --avoid AATAAA"
    )
    if not maxent_available():
        print(
            "\nNote: MaxEntScan splice scoring is unavailable; "
            "install maxentpy for --splice-backend maxent."
        )
    return 0


def _stderr(message: str) -> None:
    print(message, file=sys.stderr)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command line interface and return a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (_UsageError, SequenceError, ValueError, KeyError, OSError) as exc:
        message = exc.args[0] if exc.args else exc
        print(f"codonopt: error: {message}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
