# Legacy scripts

The two work-in-progress scripts this project started from, kept for reference
and for the benchmark comparison in `benchmarks/compare_legacy.py`. They are
not imported by the package and are not maintained.

- `optimize_codons_inprogress.py` — GC content, repeat and hairpin passes.
- `optimize_codons_with_splice_checking.py` — the same, with the hairpin pass
  replaced by MaxEntScan splice site reduction. (Previously the extension-less
  file named `in_progress_with_spliceSiteChecking`.)

Neither script contained all four passes, and both ended in a hard-coded driver
that ran the pipeline on a test protein at import time. The defects they carried
and the measured differences are listed in the main [README](../README.md#relation-to-the-original-scripts).
