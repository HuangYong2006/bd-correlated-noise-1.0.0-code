# Historical material

The current tree carries the corrected manuscript pipeline and reusable analysis/validation tools. Obsolete one-shot editing scripts, broken interfaces, superseded figure renderers, historical analysis outputs and duplicate sources were removed on 30 October 2026.

See [CLEANUP.md](CLEANUP.md) for the exact file list. All deleted tracked files remain recoverable from commit `5312322f2187bf0200e04b430f202817a94187dd` or its ancestors. Historical spectral results must not be relabelled as corrected data.

The `analyze.py`, `features.py`, `diagnostics.py` and `feature_ablation.py` tools remain because they support regression tests and reusable independent-run diagnostics. They do not establish a physical onset in the current manuscript. Optional calibration, validation and PDF-inspection utilities are retained as reusable tools rather than one-shot edits.
