"""Independent validation of required criteria, in construct coordinates."""

from __future__ import annotations

from itertools import islice

from .config import DesignConfig
from .metrics import gc_window_report, homopolymer_runs, find_hairpins
from .sequence import reverse_complement, split_codons, validate_dna
from .splice import find_splice_sites


def constraint_violations(sequence: str, config: DesignConfig) -> list[dict]:
    """Validate a CDS and return unresolved criteria with bounded region examples.

    Counts are complete. Region examples are capped at 64 per criterion; absence
    of a violation here means the configured requirement was checked in full.
    Fixed flanks are included, so a failure confined to them remains unresolved.
    """
    sequence = validate_dna(sequence, in_frame=True)
    context = config.left_flank + sequence + config.right_flank
    required = set(config.required)
    result: list[dict] = []

    def add(name: str, count: int, regions) -> None:
        if count:
            examples = list(islice(regions, 64))
            result.append({"criterion": name, "count": count, "regions": examples})

    if "motif" in required:
        regions = [
            {"start": p, "end": p + m.length, "name": m.name}
            for m in config.motifs for p in m.find_all(context)
        ]
        add("motif", len(regions), regions)
    if "gc_window" in required:
        gc = gc_window_report(context, window=config.gc_window, step=config.gc_step,
                              target_range=config.gc_range)
        add("gc_window", gc.violations,
            ({"start": p, "end": min(p + config.gc_window, len(context))}
             for p in gc.worst_positions))
    if "homopolymer" in required:
        runs = homopolymer_runs(context, config.max_homopolymer + 1)
        add("homopolymer", sum(r.length - config.max_homopolymer for r in runs),
            ({"start": r.start, "end": r.start + r.length} for r in runs))
    if "repeat" in required:
        seen: dict[str, int] = {}
        count = 0
        regions = []
        for p in range(len(context) - config.repeat_k + 1):
            kmer = context[p:p + config.repeat_k]
            if kmer in seen:
                count += 1
                if len(regions) < 64:
                    regions.extend(({"start": seen[kmer], "end": seen[kmer] + config.repeat_k},
                                    {"start": p, "end": p + config.repeat_k}))
            else:
                seen[kmer] = p
        add("repeat", count, regions)
    if "hairpin" in required:
        # Count the exact k-mer pairs used by the objective, independently.
        k = config.hairpin_stem
        count = 0
        for p in range(len(context) - k + 1):
            partner = reverse_complement(context[p:p + k])
            for q in range(p + k + config.hairpin_min_loop,
                           min(p + k + config.hairpin_max_loop + 1, len(context) - k + 1)):
                count += context[q:q + k] == partner
        examples = find_hairpins(context, stem=k, min_loop=config.hairpin_min_loop,
                                 max_loop=config.hairpin_max_loop, limit=32)
        add("hairpin", count,
            ({"start": h.stem_start, "end": h.partner_start + h.stem_length}
             for h in examples))
    if "splice" in required:
        sites = find_splice_sites(context, config.splice_scorer,
                                 donor_cutoff=config.effective_donor_cutoff,
                                 acceptor_cutoff=config.effective_acceptor_cutoff,
                                 strand=config.splice_strand)
        add("splice", len(sites),
            ({"start": s.position, "end": s.position + (9 if s.site_type == "donor" else 23),
              "type": s.site_type, "strand": s.strand} for s in sites))
    if config.excluded_codons:
        start = len(config.left_flank)
        regions = [{"start": start + 3*i, "end": start + 3*i + 3, "codon": c}
                   for i, c in enumerate(split_codons(sequence)) if c in config.excluded_codons]
        add("excluded_codons", len(regions), regions)
    return result
