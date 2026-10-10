<h1 align="center">codonopt</h1>

<p align="center">
  <strong>Codon optimization for gene synthesis and expression.</strong><br>
  Rewrite the DNA of an open reading frame so it is easier to make and express,<br>
  without changing the protein it encodes.
</p>

<p align="center">
  <a href="https://github.com/thomdehoog/computational_codon_optimization/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/thomdehoog/computational_codon_optimization/actions/workflows/ci.yml/badge.svg"></a>
  <img alt="Python 3.10+" src="https://img.shields.io/badge/python-3.10%2B-3776ab?logo=python&logoColor=white">
  <img alt="Runtime dependencies: none" src="https://img.shields.io/badge/runtime%20dependencies-none-2ea44f">
  <img alt="Interface: TypeScript" src="https://img.shields.io/badge/interface-TypeScript-3178c6?logo=typescript&logoColor=white">
  <img alt="Tests: pytest and Playwright" src="https://img.shields.io/badge/tests-pytest%20%2B%20Playwright-0a9edc?logo=pytest&logoColor=white">
  <a href="LICENSE"><img alt="License: GPL v3" src="https://img.shields.io/badge/license-GPLv3-blue"></a>
</p>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/interface-dark.png">
    <img alt="The codonopt interface after optimizing a 120-residue protein: 16 issues down to 0, with the changed codons highlighted" src="docs/interface-light.png" width="920">
  </picture>
</p>

Give it a protein, or an existing coding sequence, and it returns a synonymous
sequence that balances GC content, repeats, hairpins, homopolymers, forbidden
motifs, cryptic splice sites and host codon usage in a single score. It is a
pure-Python package with no runtime dependencies, a command line, a Python API
and a local web interface.

## Highlights

- **One objective, not four passes.** Every concern is a weighted penalty in one score, so fixing one problem cannot silently undo another.
- **Fast.** Scoring is incremental: a codon swap re-scores only what sits near it. An 800-residue protein takes 19 s, against 234 s for the scripts this replaces.
- **The protein cannot change.** The only move is swapping a codon for a synonymous one, so correctness holds by construction, not by checking afterwards.
- **Three ways in.** `codonopt` on the command line, a small Python API, and a browser or desktop interface.
- **Checked against brute force.** The tests compare the incremental score with a from-scratch recomputation and with simple reference implementations.

## Quickstart

```bash
pip install -e .

codonopt gui        # the web interface (opens in your browser)
```

Or from the command line:

```bash
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

## Web interface

`codonopt gui` starts a local server and opens the interface in a native window
if [pywebview](https://pywebview.flowrl.com/) is installed (`pip install -e ".[gui]"`),
otherwise in your browser. Paste or open a protein or DNA sequence, pick a host
table and what to keep out, and press **Optimize**. You get:

- the sequence with every changed codon highlighted, and the protein above it;
- a before/after table, a GC plot along the sequence and the score breakdown;
- copy and download for the DNA, FASTA and a report; runs can be cancelled and keep the best result so far.

It is local only: the server binds to loopback, checks Host and Origin, requires
a per-run token, sends a strict content security policy and never touches your
files. Add `--no-window --no-browser` to just serve it.

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

## Command line

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

## Python API

```python
from codonopt import DesignConfig, OptimizerSettings, get_table, optimize, resolve_motifs

config = DesignConfig(table=get_table("human"), motifs=tuple(resolve_motifs(["golden_gate"])))
result = optimize("MKWVTFISLLLLFSSAYS", config, OptimizerSettings(seed=1))
result.sequence, result.score, result.breakdown.terms
```

`measure(sequence, config)` returns what `report` prints; `optimize_dna()`
improves an existing sequence.

## How it works

Each penalty term is tied to a position (a GC window, a homopolymer, a motif
match, a splice window, a codon), so a codon swap only re-scores what sits near
it. Repeats and hairpins are not local, so an incrementally updated k-mer index
carries them. A candidate swap costs time proportional to the window sizes, not
the sequence length.

The search sweeps every codon greedily until nothing improves, optionally
followed by simulated annealing and several restarts.
`SequenceState.audit()` recomputes a score from scratch; the tests check the
incremental path against it and against brute-force implementations in
`codonopt.metrics`.

## Development

The Python package needs nothing but Python. The interface is written in
TypeScript in [`web/`](web/) and compiled with esbuild to plain HTML, CSS and
JavaScript in `src/codonopt/web/`. That output is committed, so running the
interface never needs Node; Node is only for working on it.

```bash
pip install -e ".[dev]" && pytest          # Python side

cd web && npm install
npm run build      # regenerate src/codonopt/web
npm run check      # type-check, and fail if the committed build is stale
npm test           # Playwright end-to-end tests against the real Python server
```

```
src/codonopt/   optimizer, scoring, CLI and the local server (stdlib only)
web/            TypeScript and CSS sources, build script, Playwright tests
benchmarks/     reproducible performance comparison
legacy/         the two original scripts this replaces, kept for reference
```

## History

This replaces two work-in-progress scripts, kept in [`legacy/`](legacy/). They
ran four independent passes that could undo each other, silently ignored the
local GC target (a fraction compared against a percentage range), and reloaded
the MaxEntScan matrices on every call. On an 800-residue protein the rewrite
takes 19 s instead of 234 s and ends with a better codon adaptation index
(0.90 vs 0.72) and no forbidden motifs (vs four). Details and a reproducible
benchmark are in [`legacy/README.md`](legacy/README.md).

## Limitations

- Hairpins are inverted repeats, not free-energy folds.
- The bundled `human` and `ecoli` tables are approximate; measure your own host
  for production work. The consensus splice screen is a filter, not a model.
- The search is local. A given `--seed` is reproducible, but another seed may do
  slightly better; use `--restarts` or `--anneal` when it counts.

## License

[GPL-3.0](LICENSE) © Thom de Hoog
