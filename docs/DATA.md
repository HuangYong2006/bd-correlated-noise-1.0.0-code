# Data

Simulation output is stored as compressed NumPy archives, one file per `(theta, L)` cell.
This document defines the file naming, the archive schema, the campaign inventory and how the
data are versioned.

---

## 1. Location and naming

```
data/raw/<tag>/<rule>_<theta-token>_L<L>.npz
```

The `theta` token encodes sign and decimal point so that the name is filesystem- and
shell-safe:

| `theta` | token |
|---|---|
| `0.0` | `p0d0000` |
| `1/6 = 0.1667` | `p0d1667` |
| `-0.1` | `m0d1000` |
| `0.05` | `p0d0500` |

`+` becomes `p`, `-` becomes `m`, `.` becomes `d`, and the value is formatted with four
decimals. Consequently `theta = 0.2` and `theta = 0.2000` name the same file: **the effective
precision of the `theta` axis on disk is `1e-4`.** Values closer together than that cannot
coexist in one tag.

Example: `data/raw/main/bd_p0d1667_L512.npz` is the `theta = 1/6`, `L = 512` cell of the
`main` campaign.

---

## 2. Archive schema

Every archive is written by `bd_core.save_result` and read by `bd_core.load_result`. All
files in a tag share this schema.

| Key | Type / shape | Meaning |
|---|---|---|
| `spec` | JSON string | The full `SimSpec` that produced the cell: `theta`, `L`, `T`, `n_runs`, `rule`, `n_check`, `n_w`, `seed`, `check_frac_min`, `N_ffgn`, `spectrum`. **This is the authoritative record of the cell's parameters.** |
| `N_ffgn` | `int64` scalar | Number of FFGN AR components used (`0` for the spectral generators). |
| `check_times` | `(n_check,) int64` | Sweep numbers at which full height profiles were stored. |
| `w_times` | `(n_w,) int64` | Sweep numbers of the interface-width time series. |
| `w_series` | `(n_runs, n_w) float32` | Interface width `W^2` per run, sampled at `w_times` with a fixed stride. |
| `profiles` | `(n_runs, n_check, L) int32` | Height profiles `h(x)` at each checkpoint. This is the field the diagnostics consume. |
| `stats` | `(n_runs, 3) float64` | Per-run `[mean eta', number of checkpoints captured, final W]`. For the spectral generators this array is a placeholder of zeros. |
| `wall_time` | `float64` scalar | Wall-clock seconds for the cell. |

Notes for anyone writing their own reader:

- `profiles` is the large field. For `L = 512`, `n_runs = 96` and `n_check = 32`, it is
  `96 x 32 x 512 x 4 B ≈ 6.3 MB` uncompressed per cell.
- `spec` is stored as a JSON **string**, not as a nested array, so `np.load` on a modern
  NumPy needs `allow_pickle=False` — which is what `load_result` uses. There is no pickled
  content and no arbitrary-code-execution risk in these archives.
- `check_times` is geometric: `np.geomspace(t0, T, n_check)` with
  `t0 = max(16, round(T * check_frac_min))`. Profile `[..., -1, :]` is therefore the final,
  saturated configuration, and `[..., k0:, :]` is the late-time window used by
  `analyze.py --checkpoints`.
- `w_series` stores `W^2`, not `W`. `analyze.py` takes the square root after averaging over
  runs.

Minimal reader:

```python
import sys; sys.path.insert(0, "code")
from bd_core import load_result

res = load_result("data/raw/main/bd_p0d1667_L512.npz")
print(res["spec"])                      # dict, the exact SimSpec used
print(res["profiles"].shape)            # (n_runs, n_check, L)
w = res["w_series"].mean(axis=0) ** 0.5 # ensemble-averaged interface width
```

---

## 3. Campaign inventory

| tag | generator | `spectrum` | `L` | `theta` | cells | runs / cell |
|---|---|---|---|---|---|---|
| `main` | FFGN | — | 64, 128, 256, 512 | 22 values, `0 … 0.40` | 88 | 256 |
| `spec` | spectral | `powerlaw` | 32, 64, 128, 256, 512 | 26 values, `-0.2 … 0.40` | 115 | 96 |
| `fgnspec` | spectral | `fgn` | 32, 64, 128, 256, 512 | 8 values, `0 … 0.40` | 40 | 40 |
| `explore` | mixed | mixed | 256, 1024 | 9 values | 18 | scouting |
| `testspec` | spectral | `powerlaw` | 64, 256 | 4 values | 8 | scouting |

The `spec` cell count of 115 is taken from the archived analysis header, not from the
manifest (see the caveat in section 5). `explore` and `testspec` are preliminary scouting
campaigns: they were used to choose the production grids and are not the source of any
reported number.

Derived analysis output lives in `data/analysis/<tag>_bd/`:

| File | Contents |
|---|---|
| `diagnostics.csv` | One row per `(theta, L)`: both diagnostics (`dH`, `zH`, terrace and control variants), `iForest` outputs, `alpha_s`, local exponents, morphology descriptors, `W_sat`, and the three `beta` growth windows. |
| `summary.json` | Onsets `theta*(L)` for every criterion, both finite-size-scaling fits (free and forced-zero), and the direct exponents per `theta`. |
| `rows.npz` | The same table as float arrays, for plotting. |

Plus standalone calibration/validation arrays at the top level of `data/`
(`generator_calibration*.npz`, `generator_amplitude.npz`, `generator_limit.npz`) and the
plain-text logs of analysis runs (`analysis_*.txt`), which are kept because they record the
exact console output — including the cell and grid inventory — of the published analysis.

---

## 4. Versioning policy

- **Raw `.npz` cells are not tracked by git** (see `.gitignore`). They are hundreds of MB and
  fully regenerable from `run_campaign.py` given the seeds recorded in each cell's `spec`.
- **Derived analysis tables are tracked**: `data/analysis/**/diagnostics.csv`,
  `summary.json` and `rows.npz`. They are small and they are what the manuscript numbers and
  figures are computed from, so they must be reviewable in a diff.
- **Figures are tracked** in `figs/ctp/` as PDF/SVG/TIFF together with the alignment and
  collision audit JSON, so a figure change is visible and auditable in the repository.
- For archival release, publish the full `data/raw/` tree as a **Zenodo** deposit and cite its
  DOI in the manuscript's Data Availability statement. `pack_repo.ps1 -Full` builds exactly
  this archive. Be aware of the size: the raw campaigns are **~2.5 GB uncompressed** (three
  sizes of ensemble with 32 stored profiles per run), so the full archive is ~2.5 GB and is
  suitable for Zenodo but not for a plain git push. `pack_repo.ps1` without `-Full` produces
  a 25 MB code-and-tables archive instead, which is what belongs in the repository.

---

## 5. Caveats and known issues

- **The manifests are not complete inventories.** `run_campaign.py` rewrites
  `manifest_bd.json` on every invocation, so a campaign assembled over several resumed runs
  leaves a manifest covering only the final invocation. The `spec` and `fgnspec` manifests in
  this repository are in exactly that state — the `fgnspec` manifest records a single cell.
  Use the cell files, or the `analyze.py` header line, for the real grid.
- **`stats` is zeros for the spectral generators.** `simulate_spectral` does not populate the
  per-run statistics array, so `mean eta'` cannot be read from a spectral cell. Use the FFGN
  cells or recompute it from `profiles` if you need it.
- **Cell filenames are decimal-`theta` only.** Analysis scripts recover `theta` from the
  filename by `round(int(token[3:]) / 10000.0, 4)`, which is why `spec` shows `0.1667` and
  never exactly `1/6`. Code that compares `theta` should therefore compare with a tolerance
  of `1e-4`, not for equality, whenever the value passes through a filename.
- **A killed campaign leaves truncated archives** that fail `np.load`. `verify_cells.py`
  detects them and `--delete-bad` removes them for regeneration.
