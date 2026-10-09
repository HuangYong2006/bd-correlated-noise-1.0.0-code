# Reproducing the results

Every number and figure in the manuscript is produced by the pipeline below. Commands are
written for PowerShell on Windows; the same commands work in `bash` with `python` instead of
`& $py`.

```powershell
$py = 'D:\TOOLS\anaconda3\python.exe'   # or just: $py = 'python'
```

Throughout, `--tag <name>` selects a data directory: raw cells go to `data/raw/<name>/` and
analysis output to `data/analysis/<name>_bd/`.

---

## 0. Environment

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

The first run of any script that touches `bd_core` pays a one-off numba JIT cost of roughly
20–60 s; the compiled kernels are cached in `code/__pycache__/` and reused afterwards.

Hardware used for the published results: a single workstation, 28 threads. Timings quoted
below are for that machine.

---

## 1. Validation — do this before anything else (minutes)

These are cheap and they gate every downstream interpretation.

```powershell
& $py code\check_single_step.py            # |Delta h| <= 1 and skew ≡ 0 in every cell
& $py code\check_generator_amplitude.py    # A_FFGN ~ 1.27 vs A_fGn = 2 theta (theta + 1/2)
& $py code\check_generator_limit.py        # FFGN stays correlated as theta -> 0
& $py code\calibrate_binarised.py          # theta_eff of both generators (Table: the theta_eff map)
& $py code\validate_fgn_spectrum.py        # exact-fGn filter: theta axis good to ~1%
& $py code\validate_model.py               # alpha = 1/2, beta = 1/3 at theta = 0
& $py code\compare_observables.py --tag spec --L 128   # which observable is graded in theta
```

Expected outcomes, in order: `skew` zero to machine precision; a fitted FFGN amplitude of
≈1.356 that does not move as `N` grows from 23 to 48; a `C(tau) ∝ tau^-1` tail at nominal
`theta = 0`; the `theta_eff` map of the manuscript; `alpha ≈ 0.50` and `beta ≈ 1/3` at
`theta = 0`; and the rejection of three of the four candidate observables.

`code/validate_model.py` also prints the saturation time `t_x(L)`. The campaign length is set
from it: `T = t_mult * L^1.5` with `t_mult = 15` (about 25x the measured `t_x ≈ 0.6 L^1.5`).

---

## 2. Campaigns

Three campaigns are used in the paper. Together they are roughly 5 h on 28 cores.

### 2.1 `main` — the literature-comparable FFGN axis

```powershell
& $py code\run_campaign.py --tag main --rule bd --generator ffgn `
    --L 64 128 256 512 --runs 256 --t-mult 15 `
    --thetas 0 0.0025 0.005 0.0075 0.01 0.015 0.02 0.03 0.04 0.05 0.065 0.08 `
             0.10 0.125 0.15 0.1667 0.20 0.225 0.25 0.30 0.35 0.40
```

22 exponents x 4 sizes, 256 runs per cell.

### 2.2 `spec` — the corrected (theta-continuous) axis

```powershell
& $py code\run_campaign.py --tag spec --rule bd --generator spectral --spectrum powerlaw `
    --L 32 64 128 256 512 --runs 96 --t-mult 15 `
    --thetas -0.2 -0.1 -0.05 0 0.02 0.04 0.06 0.08 0.10 0.11 0.12 0.13 0.14 0.15 `
             0.1583 0.1667 0.175 0.1833 0.1917 0.20 0.21 0.225 0.25 0.30 0.35 0.40
```

26 exponents x 5 sizes. The grid is deliberately refined around `theta_th = 1/6` so the
onset region is resolved, and it includes the three negative-`theta` controls
(`-0.2, -0.1, -0.05`) that separate "the noise changed" from "the universality class
changed".

The archived analysis of this tag records **115 cells**: not every `(theta, L)` combination
was run, and the negative-`theta` controls do not extend to every size. The analysis header
line printed by `analyze.py` lists exactly which `theta` and `L` values were found:

```
tag=spec rule=bd  115 cells
  thetas=[-0.2, -0.1, -0.05, 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.11, 0.12, 0.13, 0.14,
          0.15, 0.1583, 0.1667, 0.175, 0.1833, 0.1917, 0.2, 0.21, 0.225, 0.25, 0.3,
          0.35, 0.4]
  Ls=[32, 64, 128, 256, 512]
```

This is the string to compare against when you reproduce the campaign.

### 2.3 `fgnspec` — the exact-fGn-spectrum closure

```powershell
& $py code\run_campaign.py --tag fgnspec --rule bd --generator spectral --spectrum fgn `
    --L 32 64 128 256 512 --runs 40 `
    --thetas 0 0.10 0.1667 0.20 0.25 0.30 0.35 0.40
```

This campaign exists to close the one systematic that could have faked a flat `alpha`: the
bare power-law filter under-correlates at large `theta`. It uses a modest run count because
its purpose is to compare two filters at matched conditions, not to tighten error bars.

Each campaign is resumable: add `--skip-existing` to fill in only missing cells. Note that
`run_campaign.py` rewrites `manifest_bd.json` on **every** invocation, so a campaign built up
over several resumed runs leaves a manifest describing only the last invocation. The cell
files themselves are the source of truth; the manifest is an aid, not an inventory.

If a job was killed mid-write it leaves a truncated archive, which `verify_cells.py` will
find:

```powershell
& $py code\verify_cells.py --tag main --rule bd --delete-bad
```

---

## 3. Analysis

```powershell
& $py code\analyze.py --tag spec --rule bd --checkpoints 4 --if-runs 96 --hell-runs 96
& $py code\analyze.py --tag main --rule bd --checkpoints 4 --if-runs 96 --hell-runs 96
& $py code\analyze.py --tag fgnspec --rule bd --checkpoints 4 --if-runs 40 --hell-runs 40
```

`--if-runs` and `--hell-runs` are not cosmetic. Both diagnostics have a null floor that
depends on sample size, so if the number of runs differs between system sizes the extracted
onset drifts with statistical power instead of physics. Matching them at a fixed number
across every cell is what makes the `theta*(L)` comparison meaningful.

`--skip-iforest` skips the Isolation-Forest step for fast iteration on the distributional
diagnostic, the onsets and the exponents.

Outputs per tag: `diagnostics.csv` (one row per `(theta, L)` with all diagnostics,
exponents, `beta` windows and morphology descriptors), `summary.json` (onsets, both
finite-size-scaling fits, direct exponents) and `rows.npz` (the same table as arrays).

Then the two derived analyses:

```powershell
& $py code\fixed_theta_scaling.py --tag spec        # d(dH)/dlnL vs theta at fixed theta
& $py code\alpha_drift.py --pairs main:FFGN spec:spectral fgnspec:fGn-spectrum
```

`fixed_theta_scaling.py` is the decisive test and reads
`data/analysis/<tag>_bd/diagnostics.csv`, so `analyze.py` must have been run first.

Feature ablation (Section on the two diagnostics):

```powershell
& $py code\feature_ablation.py --tag fgnspec --rule bd --if-runs 40
```

Writes `data/analysis/ablation_fgnspec.{csv,json}`.

---

## 4. Figures

```powershell
& $py code\figures_ctp.py
```

Writes `figs/ctp/fig{1,2,3,4}*.pdf`, `.svg` and `.tiff`. Preview rasters go to
`.figpreview/` (override with the `BD_FIG_PREVIEW` environment variable) and are not part of
the submission set.

Before rebuilding the manuscript package, each figure must pass the rendered-geometry
collision audit with `verdict: PASS` and zero blocking findings. The audit helpers are an
optional external dependency; if they are not on `sys.path` the gate prints a skip message
and the export continues. See `docs/CODE_MAP.md` section 7 for `DSH_FIGURE_SKILL`.

---

## 5. Manuscript

```powershell
& $py code\check_tex.py manuscript\main.tex      # environments, refs, cites, floats, encoding
& $py code\table_width.py manuscript\main.tex    # tabular widths vs the 246 pt PRE column
& $py code\build_overleaf.py                     # rebuild manuscript/overleaf/ + upload zip
```

`pdflatex` is not required for these gates, and the manuscript cannot be compiled in an
environment without a TeX distribution. The two static checkers plus the figure collision
audit are what was used in place of a compile.

To compile the paper properly, upload `manuscript/overleaf_bd_paper.zip` to Overleaf (it runs
`pdflatex` + `bibtex` automatically), or run locally:

```bash
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

---

## 6. Determinism and repro-ducibility notes

- Seeds derive from `--seed` through `numpy.random.SeedSequence`, and the simulator's
  splitmix64 RNG is counter-based and indexed by `(run, t, x)`. Results therefore do **not**
  depend on the thread count or on the order runs are scheduled.
- The archived `manifest_bd.json` in each `data/raw/<tag>/` records, per cell, the exact
  `theta`, `L`, `T`, run count, filename and wall time, so a regenerated campaign can be
  checked cell by cell against the published one.
- Complete expected runtimes, per stage, are listed in `docs/DATA.md`.

---

## 7. Minimal end-to-end smoke test

To confirm a checkout works without committing to the full campaign:

```powershell
& $py code\validate_model.py
& $py code\run_campaign.py --tag demo --rule bd --L 64 --thetas 0 0.1667 0.40 --runs 64 --t-mult 15
& $py code\analyze.py --tag demo --rule bd --checkpoints 4 --if-runs 32 --hell-runs 32
& $py code\check_single_step.py --tag demo
```

This exercises every code path including the pipeline I/O. Onset extraction needs at least
three system sizes, so `summary.json` for `demo` will contain `null`/`NaN` fit entries — that
is expected, not a failure.
