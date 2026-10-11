"""Before/after reporting for an optimization run."""

from __future__ import annotations

from dataclasses import asdict, fields
import hashlib
import json
import platform
from importlib.metadata import version, PackageNotFoundError

from .config import DesignConfig
from .constraints import unresolved_is_all_fixed_context
from .metrics import (
    count_hairpin_pairs,
    find_hairpins,
    gc_window_report,
    homopolymer_excess,
    homopolymer_runs,
    motif_hits,
    repeat_report,
)
from .optimize import OptimizationResult
from .sequence import gc_fraction, validate_dna
from .splice import find_splice_sites

__all__ = ["measure", "report_dict", "format_report"]


def measure(sequence: str, config: DesignConfig) -> dict:
    """Measure every constrained property of ``sequence``.

    Returns a JSON-serialisable dictionary; the same keys appear for the
    starting and the optimized sequence so the two can be compared directly.
    """
    sequence = validate_dna(sequence, in_frame=True)
    context = config.left_flank + sequence + config.right_flank
    gc = gc_window_report(
        context,
        window=config.gc_window,
        step=config.gc_step,
        target_range=config.gc_range,
    )
    repeats = repeat_report(context, config.repeat_k)
    hairpins = find_hairpins(
        context,
        stem=config.hairpin_stem,
        min_loop=config.hairpin_min_loop,
        max_loop=config.hairpin_max_loop,
        limit=10,
    )
    all_runs = homopolymer_runs(context, 1)
    runs = [run for run in all_runs if run.length > config.max_homopolymer]
    hits = motif_hits(context, list(config.motifs))

    summary = {
        "length": len(sequence),
        "context_length": len(context),
        "cds_start": len(config.left_flank),
        "gc_fraction": round(gc_fraction(sequence), 4),
        "gc_windows_scored": gc.windows,
        "gc_windows_below_range": gc.low,
        "gc_windows_above_range": gc.high,
        "gc_window_min": round(gc.min_fraction, 4),
        "gc_window_max": round(gc.max_fraction, 4),
        "cai": round(config.table.cai(sequence), 4),
        "repeat_k": config.repeat_k,
        "redundant_kmers": repeats.redundant,
        "distinct_repeated_kmers": repeats.distinct_repeated,
        "repeat_covered_fraction": round(repeats.covered_fraction, 4),
        "hairpin_pairs": count_hairpin_pairs(
            context,
            stem=config.hairpin_stem,
            min_loop=config.hairpin_min_loop,
            max_loop=config.hairpin_max_loop,
        ),
        "longest_hairpin_stem": hairpins[0].stem_length if hairpins else 0,
        "homopolymer_runs_over_limit": len(runs),
        "homopolymer_excess_bases": homopolymer_excess(
            context, config.max_homopolymer
        ),
        "longest_homopolymer": max((run.length for run in all_runs), default=0),
        "forbidden_motif_hits": sum(len(positions) for positions in hits.values()),
        "forbidden_motifs_present": sorted(hits),
    }

    if config.splice_scorer is not None:
        sites = find_splice_sites(
            context,
            config.splice_scorer,
            donor_cutoff=config.effective_donor_cutoff,
            acceptor_cutoff=config.effective_acceptor_cutoff,
            strand=config.splice_strand,
        )
        summary["splice_backend"] = config.splice_scorer.name
        summary["splice_strand"] = config.splice_strand
        summary["splice_donor_sites"] = sum(
            1 for site in sites if site.site_type == "donor"
        )
        summary["splice_acceptor_sites"] = sum(
            1 for site in sites if site.site_type == "acceptor"
        )
        summary["splice_max_score"] = (
            round(max(site.score for site in sites), 2) if sites else None
        )
    return summary


def report_dict(result: OptimizationResult, config: DesignConfig) -> dict:
    """Return a complete, JSON-serialisable record of a run."""
    from . import __version__

    configuration = {f.name: getattr(config, f.name) for f in fields(config)
                     if f.name not in {"table", "motifs", "splice_scorer", "weights"}}
    configuration["weights"] = asdict(config.weights)
    configuration["motifs"] = [{"name": m.name, "pattern": m.pattern} for m in config.motifs]
    table = {"name": config.table.name, "description": config.table.description,
             "frequencies": dict(config.table.frequencies), "weights": dict(config.table.weights)}
    table["sha256"] = hashlib.sha256(json.dumps(table, sort_keys=True, allow_nan=False).encode()).hexdigest()
    configuration["table"] = table
    scorer = config.splice_scorer
    try:
        backend_version = (
            version("maxentpy") if scorer and scorer.name == "maxent"
            else __version__ if scorer and type(scorer).__module__ == "codonopt.splice"
            else None
        )
    except PackageNotFoundError:
        backend_version = "uninstalled source; record source revision separately"
    configuration["splice_scorer"] = {
        "name": scorer.name if scorer else "none",
        "class": f"{type(scorer).__module__}.{type(scorer).__qualname__}" if scorer else None,
        "version": backend_version,
        "donor_cutoff": config.effective_donor_cutoff,
        "acceptor_cutoff": config.effective_acceptor_cutoff,
    }
    return {
        "schema_version": 1,
        "software": {"codonopt": __version__, "python": platform.python_version()},
        "configuration": configuration,
        "valid": result.valid,
        "unresolved_constraints": result.unresolved_constraints,
        "stop_reason": result.stop_reason,
        "cancelled": result.cancelled,
        "notes": result.notes,
        "initial_sequence": result.initial_sequence,
        "winning_seed_sequence": result.winning_seed_sequence,
        "protein_length": len(result.protein),
        "protein": result.protein,
        "codon_usage_table": config.table.name,
        "score": {
            "initial": round(result.initial_score, 4),
            "final": round(result.score, 4),
            "improvement": round(result.improvement, 4),
            "initial_terms": {
                key: round(value, 4)
                for key, value in result.initial_breakdown.terms.items()
            },
            "final_terms": {
                key: round(value, 4) for key, value in result.breakdown.terms.items()
            },
        },
        "search": {
            "settings": result.search_settings,
            "seed": result.effective_seed,
            "restarts_run": result.restarts_run,
            "greedy_sweeps": result.sweeps,
            "moves_evaluated": result.moves_evaluated,
            "moves_accepted": result.moves_accepted,
            "seconds": round(result.elapsed, 3),
            "hit_time_limit": result.hit_time_limit,
            "cancelled": result.cancelled,
            "repair_evaluated": result.repair_evaluated,
        },
        "initial": measure(result.initial_sequence, config),
        "final": measure(result.sequence, config),
        "sequence": result.sequence,
        "context_sequence": config.left_flank + result.sequence + config.right_flank,
    }


_ROWS: tuple[tuple[str, str, str], ...] = (
    ("GC fraction", "gc_fraction", "{:.3f}"),
    ("GC windows below range", "gc_windows_below_range", "{}"),
    ("GC windows above range", "gc_windows_above_range", "{}"),
    ("Codon adaptation index", "cai", "{:.3f}"),
    ("Redundant k-mers", "redundant_kmers", "{}"),
    ("Repeat coverage", "repeat_covered_fraction", "{:.3f}"),
    ("Hairpin pairs", "hairpin_pairs", "{}"),
    ("Longest hairpin stem", "longest_hairpin_stem", "{}"),
    ("Longest homopolymer", "longest_homopolymer", "{}"),
    ("Forbidden motif hits", "forbidden_motif_hits", "{}"),
    ("Splice donor sites", "splice_donor_sites", "{}"),
    ("Splice acceptor sites", "splice_acceptor_sites", "{}"),
)


def format_report(result: OptimizationResult, config: DesignConfig) -> str:
    """Return a human-readable before/after summary."""
    before = measure(result.initial_sequence, config)
    after = measure(result.sequence, config)

    lines = [
        "Requirements        " + ("valid" if result.valid else "UNRESOLVED (not a valid design)"),
        f"Protein length      {len(result.protein)} residues",
        f"Coding sequence     {after['length']} bases",
        f"Codon usage table   {config.table.name}",
        f"Objective           {result.initial_score:.2f} -> {result.score:.2f}"
        f"  (improved by {result.improvement:.2f})",
        f"Search              {result.restarts_run} restart(s), "
        f"{result.sweeps} sweep(s), {result.moves_evaluated} moves evaluated, "
        f"{result.elapsed:.1f}s"
        + ("  [time limit reached]" if result.hit_time_limit else ""),
        f"Stop reason         {result.stop_reason}",
        "",
        f"{'Property':<24}{'before':>12}{'after':>12}",
        "-" * 48,
    ]
    for label, key, template in _ROWS:
        if key not in after:
            continue
        lines.append(
            f"{label:<24}"
            f"{template.format(before[key]):>12}"
            f"{template.format(after[key]):>12}"
        )

    present = after.get("forbidden_motifs_present")
    if present:
        lines += ["", "Remaining forbidden motifs: " + ", ".join(present)]
    if result.unresolved_constraints:
        lines += ["", "Unresolved requirements: " + ", ".join(
            _describe_unresolved(v) for v in result.unresolved_constraints)]
        if unresolved_is_all_fixed_context(result.unresolved_constraints):
            lines += [
                "Every occurrence lies wholly inside the fixed flanks, so no codon",
                "choice can remove it. The coding sequence itself meets the requirements.",
            ]
    return "\n".join(lines)


def _describe_unresolved(violation: dict) -> str:
    """One criterion, noting how much of it is beyond the optimizer's reach."""
    fixed = violation.get("fixed_context_count", 0)
    if not fixed:
        return f"{violation['criterion']} ({violation['count']})"
    return f"{violation['criterion']} ({violation['count']}, {fixed} in fixed flanks)"
