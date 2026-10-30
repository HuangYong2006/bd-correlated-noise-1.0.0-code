# PRE scientific correction implementation plan

Goal: implement the review's confirmed corrections, regenerate affected evidence without mixing it with the original data, and produce an honest updated Overleaf manuscript/package.

Authorization: user explicitly said to modify directly after the concrete review and priorities. Work in the existing project, preserving archived review and raw campaigns. No additional design approval or worktree setup is needed for this authorized in-place correction.

## Constraints
- Preserve the original spec/fgnspec/main raw data and the frozen review packet.
- Keep existing dependencies and interfaces unless a verified defect requires a change.
- Do not assume historical authors used this repository's implementation.
- Do not label a finite-lag fitted shape parameter as the RG temporal exponent.
- Report only freshly regenerated numbers in the revised manuscript.

## Task 1. Spectral production layout
Owner core agent. Files code/bd_core.py and tests/test_spectral_noise.py.
- [x] Write production API tests for temporal correlation and spatial independence for both spectra and multiple L, including chunk boundaries.
- [x] Observe failure on old implementation.
- [x] Minimal correction with output still shaped (run,time,site).
- [x] Verify tests and baseline suite, preserve deterministic seeds and add corrected layout metadata.

## Task 2. Independent statistical units
Owner statistics agent. Files code/diagnostics.py, code/analyze.py, code/feature_ablation.py, tests/test_diagnostics_runs.py.
- [x] Reproduce within-run train/calibration contamination and conditional split-error mismatch with meaningful tests.
- [x] Use one time-averaged feature vector per independent run for forest training/calibration/scoring, matching caller interfaces.
- [x] Estimate Hellinger estimator uncertainty by independent-run outer bootstrap, with accurate returned semantics. Keep null calibration distinct.
- [x] Verify tests and describe any onset covariance limits, without retaining an unsupported physical onset interpretation.

## Task 3. FFGN mathematical audit
Owner theory agent. Files code/check_generator_amplitude.py, tests/test_ffgn_covariance.py, docs/revision_20261026/theory_notes.md.
- [x] Replace diagonal-only mechanism explanation with exact shared-innovation covariance.
- [x] Test against closed-form AR cases and direct finite sums.
- [x] Distinguish uniform from Gaussian innovations, finite-window shape fitting from true asymptotic exponent, and local pmf from scale diagnostics.
- [x] Provide manuscript-ready derivation and source-grounded literature boundaries. Do not edit main.tex or bd_core.py.

## Task 4. New evidence and manuscript
Owner root. Benchmark corrected API, select and record a matched core grid, use fresh revision-specific raw/analysis tags, then regenerate figures and all affected manuscript quantities. Preserve previous data; no full-scale conclusion may use stale spectral output. Rewrite main.tex and Chinese abstract around observed results and correct claim scope. Update reproducibility instructions and Overleaf package.

## Task 5. Verification
Run full test suite, validate production time/space correlations, inspect new data provenance, compile revised LaTeX, render and inspect figures, verify zip/source identity and source/data links, and obtain independent final review. Final delivery must distinguish completed regeneration from any genuinely outstanding science.

