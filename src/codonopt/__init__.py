"""Codon optimization of an open reading frame for synthesis and expression.

A coding sequence can be written in many ways without changing the protein it
encodes, and some of those ways are much easier to synthesise and express than
others. This package resolves configured mandatory requirements before
improving soft sequence-quality and host-codon-preference objectives.

Typical use::

    from codonopt import DesignConfig, OptimizerSettings, optimize, get_table

    result = optimize(
        "MKWVTFISLLLLFSSAYS",
        DesignConfig(table=get_table("human")),
        OptimizerSettings(seed=1),
    )
    print(result.sequence, result.score)

The protein is invariant: the only move the optimizer can make is to swap one
codon for a synonymous one.
"""

from __future__ import annotations

__version__ = "0.3.0"

from .config import DesignConfig, OptimizerSettings, Weights
from .constraints import constraint_violations, unresolved_is_all_fixed_context
from .index import KmerIndex
from .metrics import (
    count_hairpin_pairs,
    find_hairpins,
    gc_window_report,
    homopolymer_excess,
    homopolymer_runs,
    repeat_report,
)
from .motifs import MOTIF_GROUPS, RESTRICTION_SITES, Motif, resolve_motifs
from .optimize import OptimizationResult, optimize, optimize_dna, seed_sequence
from .report import format_report, measure, report_dict
from .sequence import (
    SequenceError,
    gc_fraction,
    reverse_complement,
    translate,
    validate_dna,
    validate_protein,
)
from .splice import find_splice_sites, get_scorer, maxent_available
from .state import ScoreBreakdown, SequenceState
from .tables import CodonUsage, available_tables, get_table, load_table

__all__ = [
    "__version__",
    "CodonUsage",
    "DesignConfig",
    "KmerIndex",
    "MOTIF_GROUPS",
    "Motif",
    "OptimizationResult",
    "OptimizerSettings",
    "RESTRICTION_SITES",
    "ScoreBreakdown",
    "SequenceError",
    "SequenceState",
    "Weights",
    "available_tables",
    "constraint_violations",
    "count_hairpin_pairs",
    "find_hairpins",
    "find_splice_sites",
    "format_report",
    "gc_fraction",
    "gc_window_report",
    "get_scorer",
    "get_table",
    "homopolymer_excess",
    "homopolymer_runs",
    "load_table",
    "maxent_available",
    "measure",
    "optimize",
    "optimize_dna",
    "repeat_report",
    "report_dict",
    "resolve_motifs",
    "reverse_complement",
    "seed_sequence",
    "translate",
    "unresolved_is_all_fixed_context",
    "validate_dna",
    "validate_protein",
]
