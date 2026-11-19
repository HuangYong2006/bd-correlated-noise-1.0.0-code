# Reproducing the 2026-11-17 research revision

Run commands from the project root. The numerical producer is frozen at SHA256 `70b22c14eefd08a2d875e163cea93971ce95395b48b5d5903c3ef414058ad607`. Never mix cached records from an altered producer. All simulation cells are resumable; their metadata rejects incompatible parameters or source hashes. All new outputs are under `data/raw/controlled_20261117`.

The tested environment is Python 3.13.9, NumPy 2.3.5, SciPy 1.16.3, Numba 0.62.1 and Matplotlib 3.10.6. This machine uses `D:/TOOLS/anaconda3/python.exe`; the commands below use `python` for portability. Times New Roman supplies figure text.

## Reanalyse without repeating growth simulations

```powershell
python code/innovation_analysis.py --source data/raw/innovation_20261031_v3 --binary-source data/raw/innovation_binary_refined_20261031 --sizes 32 64 128 256 512 --out docs/research_refocus_20261117/reproduced_analysis
python docs/research_refocus_20261117/audit_existing.py
python code/controlled_analysis.py
python code/audit_controlled_campaign.py
python code/export_controlled_results.py
python code/refocus_figures.py
python code/controlled_figures.py
python -m pytest -q -p no:cacheprovider
```

The legacy audit reads the frozen stationary/binary archives and their analyses. The new analysis uses 2000 whole-run bootstrap draws, preserves duration pairing, and samples different sizes and G/U inputs independently. All intervals are pointwise and exploratory. Adding a source cell changes the deterministic bootstrap stream; always regenerate the final summaries and figures after the full campaign is present. The exporter refuses an incomplete size extension or control grid.

## Recreate the new experiments from seeds

These commands may take hours on a desktop. To regenerate rather than reuse checkpoints, choose a fresh output directory and analyse that directory with `--source`; do not delete historical inputs.

```powershell
python code/controlled_roughening.py --theta .25 --sizes 64 128 256 512 --runs 16 --workers 1 --seed 20261117 --out data/raw/controlled_20261117/baseline
python code/controlled_roughening.py --theta .25 --sizes 1024 --runs 64 --workers 2 --seed 20261117 --out data/raw/controlled_20261117/baseline
python code/controlled_roughening.py --theta .25 --sizes 128 256 --runs 32 --workers 1 --seed 20261118 --period-factor 4 --observation-factors 1 2 4 --out data/raw/controlled_20261117/period4
python code/controlled_roughening.py --theta .25 --sizes 128 256 --runs 32 --workers 1 --seed 20261119 --period-factor 8 --observation-factors 1 2 4 --out data/raw/controlled_20261117/period8
python code/controlled_roughening.py --theta .25 --sizes 128 --runs 32 --workers 1 --seed 20261120 --period-factor 4 --observation-factors 1 2 4 --n-modes 27 --out data/raw/controlled_20261117/modes_plus4
python code/controlled_roughening.py --theta .25 --sizes 256 --runs 32 --workers 1 --seed 20261120 --period-factor 4 --observation-factors 1 2 4 --n-modes 28 --out data/raw/controlled_20261117/modes_plus4
python code/calibration_sensitivity.py
```

For the actual L=1024 execution, the initial 16 records used two scheduling workers. `code/complete_controlled_campaign.py` resumed the frozen producer with four scheduling workers, retaining every completed record. It changes executor concurrency only; the numerical source, solver, seeds and chunk size remain fixed. `baseline/execution_L1024.json` records the launcher and producer hashes. The two-worker command above reproduces the same seeded numerical data. Wall-clock metadata can differ across machines.

## Protocols and independence

| Cohort | L | Runs per G/U input | P | Observations | N |
|---|---|---:|---|---|---|
| Baseline diagnostics | 64,128,256,512 | 16 | 2T | T,2T | N(2T) |
| Baseline extension | 1024 | 64 | 2T | T,2T | 27 |
| Period 4 | 128,256 | 32 | 4T | T,2T,4T | 23,24 |
| Period 8 | 128,256 | 32 | 8T | T,2T,4T | 23,24 |
| Modes +4 | 128,256 | 32 | 4T | T,2T,4T | 27,28 |

`T = smooth_ceil(round(15 L**1.5))`; for increasing baseline L, T is 7680, 21870, 61440, 174960 and 491520. G/U innovations have zero mean and unit variance, share the complex transfer function, and have independent seeds. Records are circular and DC removed; the interface uses each observed prefix once and never replays the record. The complete baseline/control grid contains 640 independent interface trajectories (320 saved pairs of independent G/U runs).

The calibration experiment uses 64 bootstrap refits plus the original spectrum, all propagated through the same 16 fresh innovation records. The 1040 propagated trajectories are not 1040 independent innovation records. Its percentile spread is conditional sensitivity, not a full physical confidence interval or a bound over every covariance within 0.005 of the target.

The L=64–1024 exponent curve uses historical stationary G/U ensembles (64 runs each) for L≤512 and the new 64-run L1024 ensemble. The new 16-run small-size diagnostic cohorts are not pooled into that curve. Historical innovations have variance 1/12, related to new unit-variance inputs by a common positive scale; zero-threshold signs are invariant under this scale. Both use the same period, N, time sampling and final-20% window. Old width storage used float32, new storage float64. Old stored spectra use `abs(FFT)**2/L**2`; multiplying by L gives the common normalization, checked by Parseval.

## Files and packages

- `bd_strict_binary_data_20261117.zip`: new numerical aggregates and per-run checkpoints, calibration fits, code/tests, exact main-figure source archives, frozen analysis summaries, source CSVs and hashes. It also carries the unchanged historical `innovation_controls_data_20261031.zip`, which preserves initial failed validation records and refinements. Extract only missing members of that nested historical package into the same project root when re-running the full legacy audit; retain existing new code and manifests instead of overwriting them.
- Historical 105-cell production results are preserved as their frozen summary/tables; the new package includes raw production cells needed for the displayed main figures. It does not claim to redistribute all 105 historical raw cells or to have rerun them. Supplement figures are preserved historical figures.
- `campaign_audit.json`: complete counts, unique seeds, hashes, snapshot invariants, Parseval and reaggregation checks. Noise deviations are descriptive; no new post-hoc acceptance gate is imposed.
- `figure_provenance.json` and `convergence_provenance.json`: source/output hashes and figure inference settings. Source CSVs accompany every new figure.

The local timing benchmark directory `baseline_fast` repeats one seed with another FFT chunk size. It is excluded from every scientific aggregate, figure and numerical release.

## Interpretation boundary

Time, finite period and modal truncation are separate control directions. Null-compatible duration intervals do not prove saturation or cutoff convergence. The full pointwise control table includes non-unit and null-compatible ratios; the same P=4T reference is shared across two contrasts. No multiplicity-adjusted mechanism claim is made. Likewise, binary covariance matching is approximate over measured lags, and one conditional calibration sensitivity experiment does not bound all remaining second-order effects.
