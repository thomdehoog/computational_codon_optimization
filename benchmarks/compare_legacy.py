"""Compare the legacy scripts with the package on matched operations.

Run it from the repository root::

    python benchmarks/compare_legacy.py --residues 200

The legacy scripts are executed from ``legacy/`` with their driver code
stripped, so only their function definitions run. Timings depend on the
machine; the ratios are the point.
"""

from __future__ import annotations

import argparse
import random
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from codonopt import (  # noqa: E402  (path set up above)
    DesignConfig,
    OptimizerSettings,
    get_scorer,
    get_table,
    measure,
    optimize,
    repeat_report,
    resolve_motifs,
)
from codonopt.metrics import find_hairpins
from codonopt.splice import find_splice_sites

LEGACY_SPLICE = REPO_ROOT / "legacy" / "optimize_codons_with_splice_checking.py"
LEGACY_HAIRPIN = REPO_ROOT / "legacy" / "optimize_codons_inprogress.py"
RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


def load_legacy(path: Path) -> dict:
    """Execute only the function definitions of a legacy script."""
    source = path.read_text()
    marker = 'print("STEP 1:'
    if marker in source:
        source = source[: source.index(marker)]
    namespace: dict = {}
    exec(compile(source, str(path), "exec"), namespace)  # noqa: S102
    return namespace


def timed(label: str, function, *args, **kwargs):
    start = time.perf_counter()
    result = function(*args, **kwargs)
    elapsed = time.perf_counter() - start
    return label, elapsed, result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--residues", type=int, default=200)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--skip-pipeline",
        action="store_true",
        help="skip the end-to-end comparison, which is the slow part",
    )
    args = parser.parse_args()
    # The legacy functions use the module-level generator, not the local RNG.
    random.seed(args.seed)

    splice_legacy = load_legacy(LEGACY_SPLICE)
    hairpin_legacy = load_legacy(LEGACY_HAIRPIN)
    usage = splice_legacy["codon_usage_full"]

    rng = random.Random(args.seed)
    protein = "M" + "".join(rng.choice(RESIDUES) for _ in range(args.residues - 1))
    sequence = splice_legacy["generate_codon_sequences"](
        protein, usage, max_alternatives=50, input_type="amino_acid"
    )[0]
    print(f"protein: {len(protein)} residues | coding sequence: {len(sequence)} bases\n")

    rows: list[tuple[str, float, float]] = []

    _, legacy_time, _ = timed(
        "hairpins", hairpin_legacy["find_hairpins"], sequence, 10, 20, 3
    )
    _, new_time, _ = timed(
        "hairpins", find_hairpins, sequence, stem=10, min_loop=3, max_loop=60
    )
    rows.append(("hairpin detection", legacy_time, new_time))

    _, legacy_time, _ = timed(
        "repeats (defaults)", splice_legacy["find_repeats"], sequence, 1, None
    )
    _, new_time, _ = timed("repeats", repeat_report, sequence, 12)
    rows.append(("repeat detection", legacy_time, new_time))

    scorer = get_scorer("auto")
    _, legacy_time, legacy_sites = timed(
        "splice scan", splice_legacy["find_splice_sites"], sequence, 3.0, 3.0, "both"
    )
    _, new_time, new_sites = timed(
        "splice scan",
        find_splice_sites,
        sequence,
        scorer,
        donor_cutoff=3.0,
        acceptor_cutoff=3.0,
    )
    rows.append((f"splice scan ({scorer.name})", legacy_time, new_time))
    print(
        f"splice sites found: legacy {len(legacy_sites)}, package {len(new_sites)}"
        " (the package scores both strands with cached matrices)\n"
    )

    if not args.skip_pipeline:
        print("running the legacy four-pass pipeline (this is the slow part)...")
        start = time.perf_counter()
        legacy_sequence = splice_legacy["optimize_gc_content"](
            sequence, usage, 30, 6, (30, 70), 3, 200
        )
        legacy_sequence = splice_legacy["reduce_repeats"](
            legacy_sequence, usage, 3, 12, 15, 200
        )
        legacy_sequence = splice_legacy["optimize_splice_sites"](
            legacy_sequence, usage, protein, 3.0, 3.0, 50, 2
        )
        legacy_total = time.perf_counter() - start

        config = DesignConfig(
            table=get_table("human"),
            splice_scorer=scorer,
            motifs=tuple(resolve_motifs(["golden_gate"])),
        )
        start = time.perf_counter()
        result = optimize(
            protein, config, OptimizerSettings(seed=args.seed, max_sweeps=6)
        )
        new_total = time.perf_counter() - start
        rows.append(("whole pipeline", legacy_total, new_total))

        print("\nresulting sequence quality:")
        legacy_quality = measure(legacy_sequence, config)
        new_quality = measure(result.sequence, config)
        keys = (
            "gc_fraction",
            "gc_windows_below_range",
            "gc_windows_above_range",
            "cai",
            "redundant_kmers",
            "hairpin_pairs",
            "longest_homopolymer",
            "forbidden_motif_hits",
            "splice_donor_sites",
            "splice_acceptor_sites",
        )
        print(f"  {'property':<26}{'legacy':>10}{'package':>10}")
        for key in keys:
            print(f"  {key:<26}{legacy_quality[key]:>10}{new_quality[key]:>10}")

    print(f"\n  {'operation':<26}{'legacy':>10}{'package':>10}{'speedup':>10}")
    for label, legacy_time, new_time in rows:
        ratio = legacy_time / new_time if new_time else float("inf")
        print(f"  {label:<26}{legacy_time:>9.2f}s{new_time:>9.2f}s{ratio:>9.0f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
