# computational_codon_optimization

Rewrites the coding sequence of an open reading frame to improve synthesis
success and expression, without changing the protein it encodes.

A protein can be written as many different DNA sequences, and some are far
easier to synthesise and express than others. `codonopt` searches that space
under one weighted objective:

| Criterion | Why | Default |
| --- | --- | --- |
| Overall GC | Extreme GC complicates synthesis and PCR | target 0.50 |
| Local GC | Individual stretches matter, not just the average | 30 bp windows within 0.30–0.70 |
| Direct repeats | Mis-assembly during synthesis, recombination in vivo | no repeated 12-mers |
| Hairpins | Stall synthesis and reverse transcription | no 10 bp inverted repeat, 3–60 nt loop |
| Homopolymers | Slippage in synthesis and sequencing | at most 6 identical bases |
| Forbidden motifs | Cloning sites, poly(A) signals, anything you name | none |
| Cryptic splice sites | Unwanted splicing of a transgene | MaxEntScan or a consensus screen |
| Codon usage | Rare codons slow translation | host codon adaptation index |

The protein is invariant by construction: the only move available to the
optimizer is swapping one codon for a synonymous one.

## Install

```bash
pip install -e .
```

MaxEntScan splice scoring is optional and not on PyPI:

```bash
pip install "git+https://github.com/kepbod/maxentpy.git"
```

Without it, `--splice-backend auto` falls back to a dependency-free consensus
screen; `none` switches splice scoring off and runs faster.

## Use

```bash
codonopt optimize --protein-fasta gene.faa --out gene.fna          # design
codonopt optimize --dna-fasta old.fna --out new.fna                # improve
codonopt optimize --protein-fasta gene.faa \
    --avoid golden_gate,common_cloning,AATAAA --json report.json   # constrain
codonopt report --dna-fasta gene.fna                               # measure only
codonopt tables     # codon usage tables
codonopt motifs     # known motifs, groups, IUPAC syntax
```

Each run prints a before/after summary; `--json` writes the same numbers plus
the per-term objective and the sequence.

```
Objective           598.20 -> 4.34  (improved by 593.87)

Property                      before       after
------------------------------------------------
GC fraction                    0.475       0.500
GC windows below range             4           0
Codon adaptation index         0.796       0.847
Redundant k-mers                  19           0
Longest homopolymer                7           5
Forbidden motif hits              11           0
Splice donor sites                 6           1
```

Worth knowing: `--seed` makes a run reproducible; `--restarts`, `--sweeps` and
`--anneal N` buy more search; `--seconds` caps wall-clock time;
`--weight NAME=VALUE` re-balances any term (`--weight gc_global=200` holds the
GC target harder, `usage=0` ignores codon preference); `--donor-cutoff` and
`--acceptor-cutoff` set the highest splice score treated as acceptable, so
anything the report lists is something the optimizer was paid to remove.

```python
from codonopt import DesignConfig, OptimizerSettings, get_table, optimize, resolve_motifs

config = DesignConfig(table=get_table("human"), motifs=tuple(resolve_motifs(["golden_gate"])))
result = optimize("MKWVTFISLLLLFSSAYS", config, OptimizerSettings(seed=1))
result.sequence, result.score, result.breakdown.terms
```

`measure(sequence, config)` returns what `codonopt report` prints, and
`optimize_dna()` improves an existing sequence. Custom codon tables load from
JSON or CSV via `load_table()` or `--table-file`; the bundled `human` and
`ecoli` tables are approximate, so measure your own host for production work.

## How it works

Every penalty term is anchored to a position — a GC window start, a base that
extends a homopolymer, a motif match, a splice window, a codon — so one codon
swap only re-scores the terms near the edit. Repeats and hairpins aren't local,
so an incrementally maintained k-mer index (`codonopt.index`) carries them. A
candidate swap therefore costs time proportional to the window sizes involved,
not to sequence length.

The search itself is a greedy sweep over every codon position, repeated until a
sweep changes nothing, optionally followed by simulated annealing and several
restarts. `SequenceState.audit()` recomputes the score from scratch; the tests
check the incremental path against it and against independent brute-force
implementations in `codonopt.metrics`.

## Relation to the original scripts

The two work-in-progress scripts this started from are kept in
[`legacy/`](legacy/). They ran four independent passes — random codons, GC,
repeats, splice sites — each free to undo the previous one, with nothing
measuring the result. Replacing them fixed, among others:

- **Local GC targets did nothing.** Window GC was a fraction (0–1) compared
  against a range in percent (30–70), so every window always tested as out of
  range. A percentage-looking range is now rejected rather than misread.
- **MaxEntScan matrices were reloaded on every call** (~50 ms), inside a loop
  over up to 1000 candidates per site. Cached, a window costs ~2 µs.
- **Hairpin detection scaled roughly cubically** — 0.9 s, 5.6 s and 34 s at
  600, 1200 and 2400 bp — and repeat detection enumerated every substring
  length up to the whole sequence.
- **Whole windows were rewritten rather than whole codons**, so the reading
  frame survived only when the step was a multiple of three.
- **Synonymous codons were picked uniformly at random**, discarding host codon
  preference entirely.

| Operation | legacy | now | speedup |
| --- | --- | --- | --- |
| Hairpin detection, 600 bp | 0.89 s | 0.0004 s | ~2300x |
| Repeat detection, as defaulted | 0.16 s | 0.0002 s | ~850x |
| Whole pipeline, 800 residues | 234 s | 19 s | ~12x |

The 800-residue run also ends with a higher codon adaptation index (0.90 vs
0.72) and no forbidden motifs (vs four), because every criterion is weighed at
once. Reproduce with `python benchmarks/compare_legacy.py --residues 200`.

## Development

```bash
pip install -e ".[dev]"
pytest
```

## Limitations

- Hairpins approximate secondary structure by inverted repeats, not free
  energy; use a folding tool if that matters.
- The bundled codon tables are approximate, and the consensus splice screen is
  a filter, not a model of splicing.
- The search is local: reproducible for a given `--seed`, but another seed may
  do slightly better. Use `--restarts` or `--anneal` when it counts.
