# Code map

Per-file catalogue of the `code/` directory, plus the conventions that govern it.

**Import namespace.** `code/` is used as a *flat* module directory: driver scripts do
`sys.path.insert(0, HERE)` and then `from bd_core import ...`. Nothing here is an installable
package, and every script derives the repository root as
`os.path.dirname(os.path.dirname(os.path.abspath(__file__)))`. **Do not move a file into a
subdirectory without updating that expression**, or it will silently read and write the
wrong tree.

Legend: **lib** = importable library module (no side effects); **cli** = re-runnable command-
line entry point; **script** = direct-run script, arguments at the top of the file.

---

## 1. Core model and pipeline (live)

These are the files that produce every number in the manuscript. All are free of machine-
specific paths.

| File | Kind | Purpose |
|---|---|---|
| `bd_core.py` | lib | The model: BD and RSOS rules, FFGN generator, two spectral generators, numba `prange` ensemble simulator, `.npz` I/O helpers. |
| `run_campaign.py` | cli | Campaign runner. Writes one `.npz` per `(theta, L)` cell plus `manifest_bd.json`. Args: `--tag --thetas --L --t-mult --t-cap --runs --rule --generator --spectrum --n-check --n-w --seed --skip-existing`. |
| `features.py` | lib | Amplitude-free dimensionless feature vectors (roughness exponents, normalised roughness shape, distribution moments, facet/step fractions, entropy). |
| `diagnostics.py` | lib | The two independent diagnostics: null-calibrated Hellinger distance on exact pmfs, and the Isolation-Forest novelty detector; plus onset extraction and onset finite-size-scaling fits. |
| `analyze.py` | cli | Analysis driver. Produces `diagnostics.csv`, `summary.json`, `rows.npz`. Args: `--tag --rule --checkpoints --bins --splits --if-boot --if-trees --if-runs --seed --r-frac --hell-runs --skip-iforest --outdir`. |
| `fixed_theta_scaling.py` | cli | The decisive test: is the deviation from the KPZ reference `L`-independent at fixed `theta`? Args: `--tag --rule`. |
| `alpha_drift.py` | cli | `alpha_eff(L)` from consecutive size pairs, for both generators, against the published values. Args: `--pairs` (e.g. `main:FFGN spec:spectral fgnspec:fGn-spectrum`). |
| `feature_ablation.py` | cli | Diagnostic correlation and feature-family ablation. Args: `--tag --rule --checkpoints --if-boot --if-trees --if-runs --seed --variants --outdir`. |
| `figures_ctp.py` | script | The four manuscript figures at exactly 178 mm, exported as PDF + SVG + 600 dpi TIFF. |

## 2. Generator calibration and validation (live)

| File | Kind | Purpose |
|---|---|---|
| `validate_model.py` | script | **Run this first.** Noise validation, `theta = 0` KPZ exponents (`alpha = 1/2`, `beta = 1/3`), saturation-time measurement and timings. |
| `calibrate_binarised.py` | script | The generator calibration that matters: binarised correlation shape-matched against the binarised exact-fGn family, giving `theta_eff(theta_nom)` for both generators. |
| `calibrate_spectral.py` | script | Calibration of the `theta` axis of the bare power-law spectral family. |
| `validate_fgn_spectrum.py` | script | Validates the exact-fGn spectral density filter against the bare power law (exponent accuracy and amplitude). |
| `check_generator_amplitude.py` | script | Quantitative test of the FFGN correlation amplitude vs the exact fGn amplitude, and its convergence in the number of components `N`. |
| `check_generator_limit.py` | script | Does FFGN converge to true white noise as `theta -> 0`? Shows that it does not. |
| `check_single_step.py` | cli | Structural verification across cells: `|Delta h| <= 1`, `skew(Delta h) ≡ 0`, single-parameter pmf. Args: `--tag --rule --max-cells`. |
| `compare_observables.py` | cli | Head-to-head test of the four candidate distributional observables on identical data. Args: `--tag --rule --L --splits --r-frac`. |
| `quick_hellinger_scan.py` | cli | Fast scan of the distributional diagnostic vs `theta`. Args: `--tag --rule --L --bins --splits --r-frac`. |
| `bench_spectral.py` | script | Timing benchmark for the spectral generator and simulator. |
| `verify_cells.py` | cli | Integrity check for campaign archives; `--delete-bad` removes truncated cells so `--skip-existing` regenerates them. Args: `--tag --rule --delete-bad`. |

`calibrate_generator.py` (Gaussian-domain calibration) is retained for the record but is
**superseded** by `calibrate_binarised.py`, because the model thresholds the noise and
thresholding compresses correlations. Quoting the Gaussian-domain number as `theta_eff` is
incorrect — see the module docstring of `calibrate_binarised.py`.

## 3. Legacy figure scripts (superseded)

`figures.py`, `fig_generator.py`, `fig_exponents_calibrated.py`, `fig_onsets.py`,
`fig_morphology.py`, `fig_drift.py` produced the earlier single-column draft figures now in
`figs/generator/` and `figs/physics/`. They are kept because they document how the physics
was first explored, but they are **not** used by the current manuscript. The manuscript
figures come only from `figures_ctp.py` and live only in `figs/ctp/`.

## 4. Manuscript tooling

| File | Kind | Purpose |
|---|---|---|
| `check_tex.py` | script | Static LaTeX pre-flight standing in for a compile: environment nesting, resolved `\ref`/`\cite`, `\includegraphics` targets via `\graphicspath`, math-mode balance, brace balance, encoding. Usage: `python code/check_tex.py manuscript/main.tex`. |
| `table_width.py` | script | Estimates every tabular's rendered width against the 246 pt single PRE column / 508 pt full width, so an overflowing table is caught without compiling. Usage: `python code/table_width.py manuscript/main.tex`. |
| `build_overleaf.py` | script | Rebuilds `manuscript/overleaf/` and the upload archive `manuscript/overleaf_bd_paper.zip`. Wipes the destination first so no stale figure can be uploaded. |

## 5. One-off scratch (historical, do not expect to run)

These performed a single anchored edit on a specific revision of the manuscript. They are
kept for provenance only. Several reference anchors or APIs that no longer exist, and
re-running them would fail or corrupt the current `.tex`.

`fix_layout.py` · `fix_layout2.py` · `fix_tables.py` · `fix_inline_tables.py` ·
`transpose_tables.py` · `shorten_header.py` · `fix_author.py` · `add_citations.py` ·
`add_fgn_results.py` · `finalize_fgn.py`

Additional dead/reference-only scripts:

| File | Status | Note |
|---|---|---|
| `stage1_checks.py` | dead | Imports `generate_columns` and `_spectral_amplitude`, both removed from `bd_core.py`. Historical first-stage scouting. |
| `debug_bd.py` | dead | Imports `_bd_ensemble`, removed from `bd_core.py`. |
| `check_spectrum.py` | dead | Imports `_spectral_amplitude`, removed from `bd_core.py`. Superseded by `validate_fgn_spectrum.py`. |
| `_tmp_v7.py` | scratch | Contains a hard-coded absolute path to a local PDF; requires PyMuPDF. |
| `_tmp_which_sample_count.py` | scratch | One-off question about the effective sample count per cell. |

## 6. Optional-dependency tools

These are reusable but need packages outside `requirements.txt`:

| File | Needs | Purpose |
|---|---|---|
| `inspect_pdf.py`, `layout_map.py`, `overfull.py` | `pymupdf` | Layout and overfull-box diagnostics on a compiled PDF; take the PDF path as argument. |
| `zotero_probe.py` | `pymupdf` (optional) | Reads a local Zotero library; argparse-driven, defaults to `~/Zotero`. |

## 7. Environment variables

| Variable | Effect | Default |
|---|---|---|
| `DSH_FIGURE_SKILL` | Directory containing the optional panel-alignment / collision-audit helpers used by `figures_ctp.py`. If absent, the audit gate reports itself as skipped. | `~/.dsh/skills/nature-figure` |
| `BD_FIG_PREVIEW` | Where `figures_ctp.py` writes preview rasters (review only, not part of the submission set). | `<repo>/.figpreview` |

## 8. Conventions

- Determinism: the simulator uses a counter-based splitmix64 RNG indexed by `(run, t, x)`,
  so results are reproducible and independent of thread count. Change `--seed` to resample.
- No process-based parallelism. numba `prange` threads and threaded `scipy.fft` only.
- FFT lengths should be 5-smooth; use `bd_core.smooth_ceil`.
- Paths are always derived from `__file__`. Do not introduce cwd-relative paths.
- Text files are UTF-8 without BOM; `code/check_tex.py` fails if `main.tex` is not pure
  ASCII or valid UTF-8.
