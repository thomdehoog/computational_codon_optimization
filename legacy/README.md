# Legacy scripts

The two work-in-progress scripts `codonopt` replaced, kept for reference and as
the baseline in `benchmarks/compare_legacy.py`. They are not imported by the
package and are not maintained.

- `optimize_codons_inprogress.py` — GC content, repeat and hairpin passes.
- `optimize_codons_with_splice_checking.py` — the same with the hairpin pass
  replaced by MaxEntScan splice-site reduction (previously the extension-less
  file `in_progress_with_spliceSiteChecking`).

Neither contained all four passes, and both ended in a driver that ran a test
protein at import time.

## What was wrong

They ran four passes in sequence — random codons, GC, repeats, splice sites —
each free to undo the previous one, with nothing measuring the result.

- **Local GC targets did nothing.** Window GC was a fraction (0–1) compared
  against a range in percent (30–70), so every window always tested as out of
  range.
- **MaxEntScan matrices were reloaded on every call** (~50 ms), inside a loop
  over up to 1000 candidates per site. Cached, a window costs ~2 µs.
- **Hairpin detection scaled roughly cubically**: 0.9 s, 5.6 s and 34 s at
  600, 1200 and 2400 bp.
- **Repeat detection enumerated every substring length**, up to the whole
  sequence by default.
- **A full splice scan ran per candidate**, though a codon swap only changes the
  windows it overlaps.
- **The repeat list went stale** while the sequence was being rewritten.
- **Whole windows were rewritten rather than whole codons**, so the reading
  frame survived only when the step was a multiple of three (a step of 7 or 10
  fails with "Stop codon found in the middle of the sequence"), and a sequence
  short enough for one window raised `IndexError`.
- **Synonymous codons were picked uniformly at random**, discarding host codon
  preference.
- `raise` followed by `break` left unreachable code in several places.

## Measured

| Operation | legacy | `codonopt` | speedup |
| --- | --- | --- | --- |
| Hairpin detection, 600 bp | 0.89 s | 0.0004 s | ~2300x |
| Repeat detection, as defaulted | 0.16 s | 0.0002 s | ~850x |
| Whole pipeline, 800 residues | 234 s | 19 s | ~12x |

The 800-residue run also ends with a higher codon adaptation index (0.90 vs
0.72) and no forbidden motifs (vs four), because every criterion is weighed at
once. Timings depend on the machine; the ratios are the point. Reproduce with:

```bash
python benchmarks/compare_legacy.py --residues 200
```
