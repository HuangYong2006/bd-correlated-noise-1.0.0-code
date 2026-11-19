# Binary parity, residual sensitivity and joint limits

The package contains the numerical results for *Higher-order temporal noise and finite-scale roughening in ballistic deposition*. The historical file layout, environment and controls are described in docs/research_refocus_20261117/REPRODUCE.md. Recompute from the package root using the same Python, NumPy, SciPy, Numba and Matplotlib versions.

```powershell
$env:OPENBLAS_NUM_THREADS='1'
python code/strict_binary_analysis.py
python code/refocus_figures.py --out figs/strict_binary_20261117 --notes docs/strict_binary_20261117
python code/strict_binary_figure.py
python -m pytest -q -p no:cacheprovider
```

For fresh numerical simulations, use a copy of the source and historical inputs with an empty data/raw/strict_binary_20261117 output directory. The parity and residual commands write a complete aggregate per size. The joint simulation records an atomic checkpoint per independent run and reuses matching checkpoints when resumed. All original simulation archives are inputs only.

```powershell
python code/strict_binary.py --sizes 64 128 256 512 --runs 128 --workers 2
python code/strict_binary_limits.py residual
python code/strict_binary_limits.py joint
```

- Parity input uses site-specific independent random phases, with 128 independent runs per input and size. The arrays contain exact early width snapshots and later run-resolved sampled widths, profiles, spectra and diagnostic moments. The input fields are regenerated from the archived entropy scheme and frozen source. The exact.json file contains rational expectations summed over all four-sweep input fields on three small rings.
- The residual sensitivity study uses common Gaussian innovation records for three spectra and 64 independent runs at each of three sizes. Widths, realized binary-correlation contrasts, theoretical full-lag binary correlations, fitted spectra, and fourth moments are archived. The 576 growth trajectories share 192 distinct innovation records.
- The joint study uses 32 independent runs per input over six size/setting combinations, with no record replay. Each aggregate records the actual N, record period P, duration list, bootstrap unit, run seeds and producer/solver hashes. The supplemental table shows the final settings; the JSON analysis stores all three.
- Widths and spectra in the strict experiment use 4000 whole-run resamples in analysis. Simultaneous bands cover specified early-width, late-width and spectral-band families. Local effective exponents and the displayed per-mode spectra have pointwise intervals. Separate historical controls retain their stated bootstrap settings.
- Figure 4 uses a separate figure seed (2026111721). It reads the verified analysis and the L=512 run archive. The source CSVs record displayed estimates and uncertainty bounds. Figure 1(d) estimates the connected moment with the pooled two-point products recomputed within each resample of the surrogate and uniform validation ensembles. Table SIII uses U minus surrogate throughout the binary rows.
- Figure 4 presents the exact two-time match. The long-range spatial spectra are available as `figs/refocus_20261117/fig4_structure.pdf`.
- Figure S1 uses the correctly indexed fGn-spectrum records in `data/analysis/revision_20261026/noise_validation.json`. The same-sample check with axes transposed to `(run, x, t)` is documented in [the production noise validation notes](../revision_20261026/noise_validation_notes.md). The source CSV reports one standard error for each measured value.
- The nested historical innovation package contains the surrogate calibration and validation records. At theta=0.40 and L=32, 64, and 512, the first upper bounds on the measured-lag residual were 0.010739, 0.005850, and 0.005310, respectively. The final bounds obtained from the validation ensembles described in Sec. VI of the supplement were 0.003706, 0.002553, and 0.002503. The L=512 spectrum and interface sample were held fixed for the additional validation. Extract missing members only when reproducing the entire historical analysis. Keep the new code and manifests when paths overlap.

## Numerical conventions

| Ensemble | Seed specification |
|---|---|
| Independent signs and parity blocks | `(2026111707, L, run, input)` |
| Paired surrogate spectral perturbations | `(2026111711, L, run)` |
| Joint time, period, and mode sequence | Seed domains `2026111715`, `2026111716`, and `2026111717`; the cell metadata include size, run, and input streams. |

The manuscript's Tables SV–SVII display each ratio and interval in separate columns. The numeric values are the same as in the analysis summary. To export the numeric rows for Table SIV, run `python code/export_controlled_results.py`.

The numerical ZIP is self-contained for the new analysis and figure regeneration, and retains the historical report and raw archives needed for the displayed historical comparisons. This data release does not include article PDFs, TeX/BibTeX sources, or Overleaf bundles.
