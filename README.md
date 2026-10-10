# computational_codon_optimization

Optimizes the codon sequence of an open reading frame to improve synthesis
success and expression, without changing the protein it encodes.

A protein can be written as many different coding sequences, and some are much
easier to synthesise and express than others. `codonopt` searches that space
under a single weighted objective covering:

| Property | Why it matters | Default |
| --- | --- | --- |
| Overall GC content | Extreme GC complicates synthesis and PCR | target 0.50 |
| Local GC content | Individual stretches matter, not just the average | 30 bp windows kept within 0.30–0.70 |
| Direct repeats | Mis-assembly during synthesis, recombination in vivo | no repeated 12-mers |
| Hairpins | Stall synthesis and reverse transcription | no 10 bp inverted repeats with a 3–60 nt loop |
| Homopolymer runs | Slippage during synthesis and sequencing | at most 6 identical bases |
| Forbidden motifs | Cloning sites, poly(A) signals, anything you name | none by default |
| Cryptic splice sites | Unwanted splicing of a transgene | MaxEntScan or a consensus screen |
| Codon usage | Rare codons slow translation | host codon adaptation index |

## Install

```bash
pip install -e .
```

Optional cryptic splice site scoring with MaxEntScan needs `maxentpy`, which is
not on PyPI:

```bash
pip install "git+https://github.com/kepbod/maxentpy.git"
```

Without it, `--splice-backend auto` falls back to a dependency-free consensus
screen, and `--splice-backend none` switches splice scoring off.

## Command line

```bash
# design a coding sequence for a protein
codonopt optimize --protein-fasta gene.faa --table human --out gene.fna

# keep a cloning strategy viable and avoid poly(A) signals
codonopt optimize --protein-fasta gene.faa \
    --avoid golden_gate,common_cloning,AATAAA --out gene.fna --json report.json

# improve an existing construct instead of starting fresh
codonopt optimize --dna-fasta existing.fna --out improved.fna

# measure a sequence without changing it
codonopt report --dna-fasta gene.fna

codonopt tables    # bundled codon usage tables
codonopt motifs    # known motifs, groups and IUPAC syntax
```

Every run prints a before/after summary:

```
Protein length      302 residues
Coding sequence     906 bases
Objective           598.20 -> 4.34  (improved by 593.87)

Property                      before       after
------------------------------------------------
GC fraction                    0.475       0.500
GC windows below range             4           0
Codon adaptation index         0.796       0.847
Redundant k-mers                  19           0
Hairpin pairs                      0           0
Longest homopolymer                7           5
Forbidden motif hits              11           0
Splice donor sites                 6           1
Splice acceptor sites             12           1
```

`--json` writes the same measurements, the per-term objective values and the
sequence, for use in a pipeline.

`--donor-cutoff` and `--acceptor-cutoff` set the highest splice score treated
as acceptable; a window scoring above one is both reported and penalised, so
anything the report lists is something the optimizer was paid to remove.

Useful knobs: `--seed` for reproducible runs, `--restarts` and `--sweeps` for
more search, `--anneal N` to add simulated annealing, `--seconds` for a
wall-clock budget, and `--weight NAME=VALUE` to re-balance any penalty term
(for example `--weight gc_global=200` to hold the overall GC target harder, or
`--weight usage=0` to ignore codon preference).

## Python API

```python
from codonopt import DesignConfig, OptimizerSettings, get_table, optimize, resolve_motifs

config = DesignConfig(
    table=get_table("human"),
    gc_target=0.50,
    motifs=tuple(resolve_motifs(["golden_gate", "AATAAA"])),
)
result = optimize("MKWVTFISLLLLFSSAYS", config, OptimizerSettings(seed=1))

print(result.sequence)
print(result.score, "from", result.initial_score)
print(result.breakdown.terms)        # per-term contributions
```

`measure(sequence, config)` returns the same measurements as `codonopt report`,
and `optimize_dna(sequence, ...)` improves an existing coding sequence.

Custom codon usage tables load from JSON (`{"TTT": 0.45, ...}`) or two-column
CSV via `load_table(path)` or `--table-file`. The bundled `human` and `ecoli`
tables are approximate relative synonymous codon frequencies, fine for ranking
codons; measure your own host for production work.

## How it works

Each penalty term is attributed to a position -- a GC window start, a base that
extends a homopolymer run, a motif match, a splice window, a codon -- so after
a single codon swap only the terms anchored near the edit need recomputing.
Repeats and hairpins are not local, so they are carried by an incrementally
maintained k-mer index (`codonopt.index`) that tracks how many k-mer
occurrences are redundant and how many reverse-complement pairs sit within the
loop range. A candidate swap therefore costs time proportional to the window
sizes involved rather than to the length of the sequence.

On top of that, the search is a greedy sweep over every codon position,
repeated until a whole sweep changes nothing, optionally followed by simulated
annealing to escape a local minimum, optionally across several restarts. The
protein is invariant by construction: the only available move swaps one codon
for a synonymous one.

`SequenceState.audit()` recomputes the score from scratch, and the test suite
checks the incremental bookkeeping against it and against the independent
brute-force measurements in `codonopt.metrics`.

## Relation to the original scripts

The two work-in-progress scripts this repository started from are kept under
[`legacy/`](legacy/) for reference. They ran four independent passes -- random
codons, GC content, repeats, splice sites -- each free to undo the previous
one, with no overall measure of the result. Replacing them fixed a number of
defects:

- **Local GC targets were ignored.** Window GC was computed as a fraction
  (0–1) and compared against a range given in percent (30–70), so every window
  always tested as out of range and the requested range had no effect.
  Fractions are now used throughout, and a range that looks like percentages is
  rejected with a message rather than silently misread.
- **MaxEntScan matrices were reloaded on every call** -- about 50 ms each time,
  inside a loop over up to 1000 candidate sequences per splice site. Scoring a
  window with cached matrices takes about 2 µs rather than 6 ms.
- **Hairpin detection scaled roughly cubically** (0.9 s, 5.6 s and 34 s for
  600, 1200 and 2400 bp) by searching for a reverse complement inside a loop
  over stem sizes, start positions and loop sizes. Indexing k-mers once makes
  it linear for a fixed loop range.
- **Repeat detection enumerated every substring length** between `min_length`
  and `max_length`, building a dictionary quadratic in sequence length, with
  `max_length=None` defaulting to the whole sequence. One k-mer size in a
  single pass finds the same repeats.
- **A full splice site scan ran per candidate sequence**, although a codon swap
  can only change the windows it overlaps.
- **The repeat list went stale.** `reduce_repeats` collected repeats and their
  positions once per iteration and then kept working from that list while
  rewriting the sequence, so later replacements were chosen for content that
  had already changed.
- **Whole windows were rewritten rather than whole codons**, so `optimize_gc_content`
  only held the reading frame when the step happened to be a multiple of three.
  Nothing checked: a step of 7 or 10 fails with "Stop codon found in the middle
  of the sequence". It also raised `IndexError` on a sequence short enough to
  yield a single window.
- **`raise` followed by `break`** left unreachable code in several places, and
  `found_valid_replacement` was assigned but never read in `reduce_repeats`.
- **Synonymous codons were chosen uniformly at random**, discarding host codon
  preference; codon adaptation is now part of the objective.
- Hairpin reporting could describe one structure with a non-maximal stem,
  depending on which k-mer register the scan reached first.

Measured on the same machine, 200-residue protein unless stated:

| Operation | legacy | package | speedup |
| --- | --- | --- | --- |
| Hairpin detection | 0.89 s | 0.0004 s | ~2300x |
| Repeat detection (as defaulted) | 0.16 s | 0.0002 s | ~850x |
| Whole pipeline, 800 residues | 234 s | 19 s | ~12x |

The 800-residue comparison also ends with a higher codon adaptation index
(0.90 against 0.72) and no forbidden motifs (against four), because the
objective weighs every criterion at once instead of optimising them in
sequence.

Reproduce with:

```bash
python benchmarks/compare_legacy.py --residues 200
```

## Development

```bash
pip install -e ".[dev]"
pytest
```

The suite covers the genetic code and sequence primitives, codon usage tables
and CAI, motif resolution on both strands with IUPAC codes, the measurement
functions, the incremental index and score against brute-force recomputation
under randomised edits, protein invariance and the behaviour of each penalty
term, both splice backends, FASTA handling and the command line.

## Limitations

- mRNA secondary structure is approximated by inverted repeats, not free
  energy; for a folding-based score, pass sequences to a dedicated tool.
- The bundled codon usage tables are approximate, and the consensus splice
  screen is a filter rather than a model of splicing.
- The search is a local one. It is reproducible for a given `--seed`, but a
  different seed may find a slightly better sequence; use `--restarts` or
  `--anneal` when it matters.
