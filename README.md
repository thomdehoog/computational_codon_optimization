# computational_codon_optimization

Rewrites the DNA of an open reading frame so it is easier to synthesise and
express, without changing the protein it encodes.

Give it a protein (or an existing coding sequence) and it returns a synonymous
sequence that balances GC content, repeats, hairpins, homopolymers, forbidden
motifs, cryptic splice sites and host codon usage in a single score.

## Quickstart

```bash
pip install -e .

codonopt optimize --add-stop --avoid golden_gate,EcoRI --seed 5 --splice-backend motif \
    --protein MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEENFKALVLIAFAQYLQQCPF
```

The sequence goes to stdout as FASTA (or `--out FILE`), the summary to stderr:

```
>protein codonopt table=human score=1.17
ATGAAATGGGTGACATTCATCAGCTTACTGCTGCTGTTCTCTTCCGCTTACTCCCGGGGC
GTGTTTAGACGGGACACTCATAAATCTGAGATCGCTCACAGGTTTAAGGACCTGGGCGAG
GAGAACTTCAAGGCCTTGGTTTTAATCGCCTTCGCACAGTATCTGCAGCAGTGCCCTTTC
TGA

Objective           54.44 -> 1.17  (improved by 53.28)

Property                      before       after
------------------------------------------------
GC fraction                    0.481       0.503
Codon adaptation index         0.774       0.837
Forbidden motif hits               1           0
Splice donor sites                 2           0
Splice acceptor sites              1           0
```

`--json FILE` writes the same numbers, the per-term score and the sequence.

## What it optimizes

| Criterion | Why | Default |
| --- | --- | --- |
| Overall GC | Extreme GC complicates synthesis and PCR | target 0.50 |
| Local GC | Stretches matter, not just the average | 30 bp windows within 0.30–0.70 |
| Direct repeats | Mis-assembly in synthesis, recombination in vivo | no repeated 12-mers |
| Hairpins | Stall synthesis and reverse transcription | no 10 bp inverted repeat, 3–60 nt loop |
| Homopolymers | Slippage in synthesis and sequencing | at most 6 identical bases |
| Forbidden motifs | Cloning sites, poly(A) signals, anything you name | none |
| Cryptic splice sites | Unwanted splicing of a transgene | MaxEntScan or a consensus screen |
| Codon usage | Rare codons slow translation | host codon adaptation index |

The protein cannot change: the only move the optimizer has is swapping a codon
for a synonymous one.

## Commands

```bash
codonopt optimize --protein-fasta gene.faa --out gene.fna   # design from a protein
codonopt optimize --dna-fasta old.fna --out new.fna         # improve a sequence
codonopt report --dna-fasta gene.fna                        # measure, change nothing
codonopt tables                                             # codon usage tables
codonopt motifs                                             # known motifs and groups
```

| Option | Effect |
| --- | --- |
| `--avoid SPEC` | Exclude motifs: a group (`golden_gate`, `common_cloning`, `regulatory`), an enzyme (`EcoRI`), an IUPAC pattern (`GGWCC`) or `NAME=PATTERN`. Both strands are checked. |
| `--table NAME` / `--table-file F` | Host codon usage: `human`, `ecoli`, `uniform`, or your own JSON/CSV. |
| `--weight NAME=VALUE` | Re-balance a term: `gc_global=200` holds the GC target harder, `usage=0` ignores codon preference. |
| `--splice-backend` | `auto` (MaxEntScan if installed, else a consensus screen), `maxent`, `motif`, `none`. |
| `--seed N` | Reproducible runs. |
| `--restarts N`, `--sweeps N`, `--anneal N` | More search effort. |
| `--seconds N` | Wall-clock cap (default 120). |

`codonopt optimize --help` lists the rest. MaxEntScan is optional and not on
PyPI: `pip install "git+https://github.com/kepbod/maxentpy.git"`.

## Python

```python
from codonopt import DesignConfig, OptimizerSettings, get_table, optimize, resolve_motifs

config = DesignConfig(table=get_table("human"), motifs=tuple(resolve_motifs(["golden_gate"])))
result = optimize("MKWVTFISLLLLFSSAYS", config, OptimizerSettings(seed=1))
result.sequence, result.score, result.breakdown.terms
```

`measure(sequence, config)` returns what `report` prints; `optimize_dna()`
improves an existing sequence.

## How it works

Each penalty term is tied to a position — a GC window, a homopolymer, a motif
match, a splice window, a codon — so a codon swap only re-scores what sits near
it. Repeats and hairpins are not local, so an incrementally updated k-mer index
carries them. A candidate swap costs time proportional to the window sizes, not
the sequence length.

The search sweeps every codon greedily until nothing improves, optionally
followed by simulated annealing and several restarts.
`SequenceState.audit()` recomputes a score from scratch; the tests check the
incremental path against it and against brute-force implementations in
`codonopt.metrics`.

## History

This replaces two work-in-progress scripts, kept in [`legacy/`](legacy/). They
ran four independent passes that could undo each other, silently ignored the
local GC target (a fraction compared against a percentage range), and reloaded
the MaxEntScan matrices on every call. On an 800-residue protein the rewrite
takes 19 s instead of 234 s and ends with a better codon adaptation index
(0.90 vs 0.72) and no forbidden motifs (vs four). Details and a reproducible
benchmark are in [`legacy/README.md`](legacy/README.md).

## Development

```bash
pip install -e ".[dev]"
pytest
```

## Limitations

- Hairpins are inverted repeats, not free-energy folds.
- The bundled `human` and `ecoli` tables are approximate; measure your own host
  for production work. The consensus splice screen is a filter, not a model.
- The search is local. A given `--seed` is reproducible, but another seed may do
  slightly better; use `--restarts` or `--anneal` when it counts.
