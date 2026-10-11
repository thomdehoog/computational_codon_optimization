"""Independent validation of required criteria, in construct coordinates."""

from __future__ import annotations

from .config import DesignConfig
from .metrics import gc_window_report, homopolymer_runs, find_hairpins, window_starts
from .sequence import reverse_complement, split_codons, validate_dna
from .splice import find_splice_sites

_MAX_EXAMPLES = 64


def constraint_violations(sequence: str, config: DesignConfig) -> list[dict]:
    """Validate a CDS and return unresolved criteria with bounded region examples.

    Counts are complete. Region examples are capped at 64 per criterion; absence
    of a violation here means the configured requirement was checked in full.

    Fixed flanks are included, so a failure confined to them remains unresolved:
    the construct really does contain it. Such a failure is also beyond the
    optimizer's reach, because no codon choice can touch it, so each criterion
    additionally reports ``fixed_context_count`` -- how much of ``count`` lies
    wholly outside the coding sequence -- and every listed region carries a
    ``fixed_context`` flag. When the two counts are equal the design is as good
    as the supplied context allows, which is a different situation from a CDS
    the search failed to clean up.
    """
    sequence = validate_dna(sequence, in_frame=True)
    context = config.left_flank + sequence + config.right_flank
    cds_start = len(config.left_flank)
    cds_end = cds_start + len(sequence)
    required = set(config.required)
    result: list[dict] = []

    def editable(region: dict) -> bool:
        """True if the region overlaps the coding sequence, so a swap can reach it."""
        return region["start"] < cds_end and region["end"] > cds_start

    def add(name: str, occurrences, examples=None) -> None:
        """Record a criterion from every occurrence, optionally with chosen examples.

        ``occurrences`` yields ``(region, weight)``; the weights sum to the
        reported count, so a criterion counted in bases (homopolymers) stays
        consistent with one counted in sites.
        """
        count = fixed = 0
        found: list[dict] = []
        for region, weight in occurrences:
            count += weight
            reachable = editable(region)
            if not reachable:
                fixed += weight
            if examples is None and len(found) < _MAX_EXAMPLES:
                found.append({**region, "fixed_context": not reachable})
        if not count:
            return
        if examples is not None:
            found = [{**region, "fixed_context": not editable(region)}
                     for region in examples][:_MAX_EXAMPLES]
        result.append({
            "criterion": name,
            "count": count,
            "fixed_context_count": fixed,
            "regions": found,
        })

    if "motif" in required:
        add("motif", (
            ({"start": p, "end": p + m.length, "name": m.name}, 1)
            for m in config.motifs for p in m.find_all(context)
        ))
    if "gc_window" in required:
        low, high = config.gc_range
        def bad_windows():
            for start in window_starts(len(context), config.gc_window, config.gc_step):
                end = min(start + config.gc_window, len(context))
                fraction = (context.count("G", start, end)
                            + context.count("C", start, end)) / (end - start)
                if fraction < low or fraction > high:
                    yield {"start": start, "end": end, "gc_fraction": round(fraction, 4)}, 1
        # The worst windows make better examples than the first ones found.
        report = gc_window_report(context, window=config.gc_window, step=config.gc_step,
                                  target_range=config.gc_range)
        add("gc_window", bad_windows(),
            examples=[{"start": p, "end": min(p + config.gc_window, len(context))}
                      for p in report.worst_positions])
    if "homopolymer" in required:
        add("homopolymer", (
            ({"start": r.start, "end": r.start + r.length, "base": r.base},
             r.length - config.max_homopolymer)
            for r in homopolymer_runs(context, config.max_homopolymer + 1)
        ))
    if "repeat" in required:
        k = config.repeat_k
        seen: dict[str, int] = {}
        repeats: list[tuple[dict, int]] = []
        examples: list[dict] = []
        for p in range(len(context) - k + 1):
            kmer = context[p:p + k]
            first = seen.get(kmer)
            if first is None:
                seen[kmer] = p
                continue
            # One redundant occurrence is one violation; show both copies, since
            # either of them may be the one worth editing.
            repeats.append(({"start": p, "end": p + k, "first_start": first}, 1))
            if len(examples) < _MAX_EXAMPLES:
                examples.extend(({"start": first, "end": first + k},
                                 {"start": p, "end": p + k}))
        add("repeat", repeats, examples=examples)
    if "hairpin" in required:
        # Count the exact k-mer pairs used by the objective, independently.
        k = config.hairpin_stem
        def pairs():
            for p in range(len(context) - k + 1):
                partner = reverse_complement(context[p:p + k])
                for q in range(p + k + config.hairpin_min_loop,
                               min(p + k + config.hairpin_max_loop + 1, len(context) - k + 1)):
                    if context[q:q + k] == partner:
                        yield {"start": p, "end": q + k}, 1
        strongest = find_hairpins(context, stem=k, min_loop=config.hairpin_min_loop,
                                  max_loop=config.hairpin_max_loop, limit=32)
        add("hairpin", pairs(),
            examples=[{"start": h.stem_start, "end": h.partner_start + h.stem_length}
                      for h in strongest])
    if "splice" in required:
        sites = find_splice_sites(context, config.splice_scorer,
                                 donor_cutoff=config.effective_donor_cutoff,
                                 acceptor_cutoff=config.effective_acceptor_cutoff,
                                 strand=config.splice_strand)
        add("splice", (
            ({"start": s.position,
              "end": s.position + (9 if s.site_type == "donor" else 23),
              "type": s.site_type, "strand": s.strand}, 1)
            for s in sites
        ))
    if config.excluded_codons:
        add("excluded_codons", (
            ({"start": cds_start + 3 * i, "end": cds_start + 3 * i + 3, "codon": c}, 1)
            for i, c in enumerate(split_codons(sequence)) if c in config.excluded_codons
        ))
    return result


def unresolved_is_all_fixed_context(violations: list[dict]) -> bool:
    """True when nothing unresolved can be reached by any codon choice.

    Every remaining occurrence lies wholly in the immutable flanks, so the
    coding sequence is already as clean as the supplied construct permits.
    """
    return bool(violations) and all(
        v["count"] == v["fixed_context_count"] for v in violations
    )
