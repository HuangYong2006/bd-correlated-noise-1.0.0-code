> As of 2026-10-30, the three corrected `rev20261026_*` production campaigns (105 NPZ cells) are versioned in Git. Other raw campaigns remain excluded. New control records are documented in `docs/submission_20261030/REPRODUCE.md`.

# Data for the scientific revision

The revised manuscript uses only the three `rev20261026_*` production tracks,
the revision analysis directory and the production-noise validation JSON. Earlier
tracks are preserved as historical records, not silently relabelled as corrected.

## Locations and naming

Raw data are compressed NumPy archives at
`data/raw/<tag>/bd_<theta-token>_L<L>.npz`. For example,
`bd_p0d2500_L256.npz` is theta=0.25, L=256. Filenames round theta to four decimals;
the JSON `spec` inside each file is authoritative, including the actual 1/6 value.
All archives load with `allow_pickle=False`.

## Current production grid

| Tag | Input | Cells | Sizes | Theta | Independent runs/cell |
|---|---|---:|---|---|---:|
| rev20261026_shared | uniform-driven shared-innovation AR, iid shortcut at zero | 35 | 32,64,128,256,512 | 0,.10,1/6,.20,.25,.30,.40 | 64 |
| rev20261026_powerlaw | corrected Gaussian power-law spectral field | 35 | same | same | 64 |
| rev20261026_fgn | corrected Gaussian fGn-spectrum field | 35 | same | same | 64 |

The 105 cells contain 6720 independent runs. T is round(15*L^1.5), rounded upward
to a 5-smooth length for spectral inputs. These exact lengths are stored in spec.
The two spectral filters share matched Gaussian base streams and must not be
called statistically independent replications. The AR path uses different
uniform innovations; its numerical path was not changed by the layout fix.

## Raw schema

| Field | Shape/type | Meaning |
|---|---|---|
| spec | JSON string | Full SimSpec, including theta, L, T, n_runs, rule, seed and spectrum where applicable |
| N_ffgn | scalar integer | Number of AR modes; zero for spectral data |
| noise_layout | scalar string | `spectral-run-time-site-v2` for corrected spectral production; old missing values load as `unversioned` |
| check_times | n_check integers | Times of the geometric profile snapshots |
| w_times | n_w integers | Width sampling times |
| w_series | n_runs by n_w, float32 | Per-run squared interface width W^2, not W |
| profiles | n_runs by n_check by L, int32 | Stored height configurations; the last is the final state, not automatically a proven stationary state |
| stats | n_runs by 3 | AR bookkeeping; spectral entries remain zero placeholders |
| wall_time | scalar float | Cell wall time in seconds |

The spectral stats field cannot be used as a measured deposition fraction, and
noise correlations cannot be recovered from final height profiles. Use the
separate final-interface noise validation instead.

## Derived current evidence

`data/analysis/revision_20261026/summary.json` contains the declared protocol,
105 cell summaries, 21 generator/theta width fits and the full raw-source hashes
and parameter records. `observables.csv` and `exponents.csv` are readable tables.
The analysis explicitly rejects unversioned spectral fields and an incomplete
or unexpected production grid.

Widths average W^2 over the final 20% of stored times per run, then take the
ensemble RMS. The preceding 20% is a duration diagnostic. Widths, p_flat and
size slopes have 500-draw independent-run bootstrap SD. Hellinger errors use
200 outer run bootstraps of a complete estimator containing 120 conditional
partitions. Null-floor means/scatters describe those partitions, not estimator
errors or calibrated p values. DeltaH(0,L)=0 is a reference definition.

`noise_validation.json` contains eight validation cells, individual-run temporal
and spatial correlations, run SE, actual finite-period targets, engineering
tolerances, source hashes and a same-sample faulty-layout control. Its AR
covariances are exact stationary raw-field sums, not binary correlations.

`figs/revision_20261026/figure_sources.json` records the figure inputs and
illustrative-profile indices. All quantitative morphology panels use all 64
runs. `manuscript/revision_results.tex` is generated from current source tables.

## Legacy boundaries

- `main` has 88 cells with 256 runs each. It is an archive of this repository's
  shared-innovation variant, not proof of historical authors' exact algorithms.
- `spec` has 115 cells and `fgnspec` 40 cells. Their stored counts are 192 runs
  for L<=256 and 96 for L=512. These old spectral tracks contain the time/site
  layout error and are excluded from current results.
- `explore` and `testspec`, old calibration arrays at data/*.npz, and derived
  main_bd/spec_bd/fgnspec_bd tables are historical. They do not feed the new
  manuscript or figures.
- The pre-revision manuscript zip and frozen review packet are preserved.
  Earlier submission-readiness notes concern the superseded draft.

Raw campaigns are intentionally ignored by Git. Derived current tables and
figure sources can be versioned. Do not invent a public data DOI or publication
status; provide an actual deposit identifier only when one exists.
