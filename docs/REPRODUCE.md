> Current 2026-10-30 supplement and repository data instructions: [submission reproduction guide](submission_20261030/REPRODUCE.md). The production commands below remain applicable.

# Reproducing the revised results

This workflow reproduces the scientific revision dated 2026-10-26. It supersedes
older commands and readiness notes. The archived `main`, `spec` and `fgnspec`
tracks are not inputs to the revised manuscript and are retained only in Git history. In particular, the old spectral
tracks contain the documented time/site layout error.

## Environment and regression checks

The tested local interpreter is Python 3.13.9 at `D:/TOOLS/anaconda3/python.exe`,
with NumPy 2.3.5, SciPy 1.16.3, Numba 0.62.1 and scikit-learn 1.7.2. On another
machine use Python >=3.10 and install requirements.txt plus pytest. No new
scientific dependency was introduced by the correction.

```powershell
$py = 'D:/TOOLS/anaconda3/python.exe'
$env:NUMBA_NUM_THREADS = '12'
$env:OPENBLAS_NUM_THREADS = '1'
& $py -m pytest tests -q --basetemp .dsh_tmp/pytest
& $py code/validate_production_noise.py
& $py code/check_generator_amplitude.py
```

The production validation uses eight independent runs at theta=0.25,0.40,
L=32,64 and T=65536 for each spectral filter. Temporal targets are computed from
the actual finite-period spectrum, with its DC and alias conventions. It also
checks the same-sample axis-permuted fault control. The output is
`data/analysis/revision_20261026/noise_validation.json`, including individual-run
correlations, standard errors, engineering check tolerances and source hashes.
The covariance script prints deterministic shared/independent AR results and
does not replace old empirical amplitude archives.

## Regenerated campaigns

Each current track has 35 cells and 64 independent runs per cell. The grid was
fixed before evaluating the corrected results: theta=0,0.10,1/6,0.20,0.25,0.30,0.40
and L=32,64,128,256,512. These values cover the white reference, small memory,
the continuum reference boundary, the commonly discussed quarter range, and
stronger correlations. They are not a fine grid for measuring an onset.

```powershell
& $py code/run_campaign.py --tag rev20261026_shared --generator ffgn --rule bd --L 32 64 128 256 512 --runs 64 --t-mult 15 --thetas 0 0.10 0.1666666666666667 0.20 0.25 0.30 0.40 --skip-existing
& $py code/run_campaign.py --tag rev20261026_powerlaw --generator spectral --spectrum powerlaw --rule bd --L 32 64 128 256 512 --runs 64 --t-mult 15 --thetas 0 0.10 0.1666666666666667 0.20 0.25 0.30 0.40 --skip-existing
& $py code/run_campaign.py --tag rev20261026_fgn --generator spectral --spectrum fgn --rule bd --L 32 64 128 256 512 --runs 64 --t-mult 15 --thetas 0 0.10 0.1666666666666667 0.20 0.25 0.30 0.40 --skip-existing
```

These commands do not overwrite the original tracks. Seeds are deterministically
assigned from the base seed, theta and L. Spectral cells must record
`noise_layout=spectral-run-time-site-v2`. Counts, sizes and nominal time factor
are matched. The spectral T is rounded upward to a 5-smooth FFT length; exact T
and all parameters are stored per cell. The shared AR path is the specified
uniform-driven common-innovation variant, with an iid shortcut at theta=0;
it is not asserted to reproduce every historical FFGN implementation.

The manifest produced by run_campaign.py describes the latest invocation, which
may contain only resumed cells. The raw cell parameters and the full source list
written by revision_analysis.py are the authoritative inventory.

## Analysis, figures and manuscript

```powershell
& $py code/revision_analysis.py --bootstrap 500 --hell-bootstrap 200 --splits 120
& $py code/export_revision_tables.py
& $py code/revision_figures.py
& $py code/submission_controls.py
& $py code/build_overleaf.py
```

The analysis requires the full expected grid and refuses unversioned spectral
data. It writes summary.json, observables.csv and exponents.csv under
`data/analysis/revision_20261026/`, with raw-file hashes and exact cell specs.
The width is the ensemble RMS after averaging W^2 within each run over the final
20% of stored times. Five-size and adjacent-size slopes use an independent-run
bootstrap with 500 draws; the adjacent preceding 20% provides a duration check.
The Hellinger excess uses 120 conditional partitions and an outer 200-draw run
bootstrap of the full estimator. Conditional null scatter is not estimator error.
DeltaH at theta=0 is zero by reference definition. No theta_c fit is inferred.

LaTeX numbers are generated into manuscript/revision_results.tex. Figures go to
`figs/revision_20261026/` as PDF, SVG and 600 dpi PNG. Their source manifest records
input hashes, all-run quantitative sources and the illustrative-profile selection.
Use revision_figures.py directly; superseded plotting entry points were removed.

Build produces manuscript/overleaf_bd_paper.zip from an explicit whitelist.
Upload that zip to Overleaf, choose main.tex, and use pdfLaTeX/BibTeX. The old
draft files are no longer distributed in the current tree. The builder writes
the ZIP directly without retaining a duplicate source directory.

## Statistical tools retained outside the main inference

The legacy analyze.py and feature_ablation.py entry points were corrected to use
one late-time mean feature vector per independent run. Forest training, threshold
calibration and a held-out null use disjoint runs. Small calibration sets have
coarse FPR resolution and zero null variance yields an unavailable z score,
not artificial precision. These secondary tools are not evidence for the revised
manuscript's physical onset claims. No such claim is made.
