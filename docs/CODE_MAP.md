# Current code map

All scripts use `code/` as a flat module directory. Run them from the repository root. The canonical manuscript pipeline is:

| File | Role |
|---|---|
| `bd_core.py` | Parallel growth rules, specified AR/spectral input generation and archive I/O |
| `run_campaign.py` | Generate production cells with recorded parameters and seeds |
| `validate_production_noise.py` | Validate final binary temporal/spatial fields and layout control |
| `check_generator_amplitude.py` | Exact shared/independent innovation covariance calculations |
| `revision_analysis.py` | Independent-run widths, finite-size slopes and local distances |
| `export_revision_tables.py` | Generate manuscript numerical macros |
| `revision_figures.py` | Four production figures and the shared Times New Roman style |
| `submission_controls.py` | Paired filter bootstrap, fixed-record duration controls and supplemental figure |
| `build_overleaf.py` | Build the source ZIP directly from canonical sources |

`features.py`, `diagnostics.py`, `analyze.py` and `feature_ablation.py` provide reusable run-level diagnostics and support regression tests. Their optional novelty/onset outputs are not physical critical-point evidence.

Other retained scripts provide calibration studies, integrity checks, benchmarks, model validation, TeX checks and optional PDF/library inspection. They are outside the main evidence pipeline. Removed one-shot/broken/historical renderers are listed in [CLEANUP.md](CLEANUP.md).

[The current reproduction guide](submission_20261030/REPRODUCE.md) gives the execution order. The source package contains `main.tex`, both numerical macro files, `refs.bib`, a reading abstract, five PDF/SVG figures and a short README. It does not require a second local `overleaf/` tree.
