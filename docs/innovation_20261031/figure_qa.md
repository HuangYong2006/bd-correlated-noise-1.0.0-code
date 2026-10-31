# Figure contract and QA

Backend: Python/matplotlib, matching the existing project workflow. Preference saved locally in figure_backend.json.
Purpose: research report figures and editable source for later manuscript selection, not a claim of final PRE submission formatting.
Archetype: quantitative comparison grids.
Primary question: does matching the actual binary input covariance remove observable differences in finite-size roughness, while raw-covariance matching alone may fail after thresholding?
Hero: slope_differences, whole-run bootstrap pointwise95% intervals at both T and2T. Shared independent comparisons, raw-matched controls, binary-matched comparison are distinct rows. Exact Gaussian equality is a control, not a novelty claim.
Support: innovation_controls shows widths with95% run-bootstrap intervals and input correlation curves; structure_factors shows final spatial spectra at the largest completed size. Noise/structure curves are descriptive ensemble means; uncertainty and matching validity are reported separately in summary/report tables.
Export: editable SVG/PDF plus PNG previews; report-sized panels retain >=5pt text. Select final manuscript physical dimensions after evidence selection.
Data integrity: all final v3 formal cells, all seven variants, both durations. v1/v2 trial archives excluded for documented implementation/provenance reasons. Partial plots are explicitly monitoring artifacts.
Caption requirements: n=64 independent interface runs per variant, held-out uniform reference n=32 for binary comparisons,20% late-time window of regularly stored samples;2000 run-bootstrap draws; pointwise intervals, exploratory comparisons; conditional on fitted binary spectrum. Structure-factor normalization sum_k S(k)=Var_x(h); no fitted spectral exponent.
Pending final QA: inspect all rendered figures, source preflight, PDF glyph audit, and verify hashes/counts after both formal jobs finish.

Visual QA of interim slope figure: readable labels, interval bars, no clipping. Missing binary comparison rows in partial monitoring reflect unfinished matching cells; final figure must have all15 controls before delivery. Largest-size input and structure-factor panels must be re-inspected at completion. Glyph floor on the interim PDF is8pt. Static preflight warnings are expected for report-size vector+PNG exports; RNG is run bootstrap, not invented observations. Mean squared widths must be checked positive before interpreting log-slope fits.
