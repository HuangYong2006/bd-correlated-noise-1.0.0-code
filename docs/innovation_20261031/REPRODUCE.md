# Reproduce the current innovation-control revision

The current manuscript uses the original 105 corrected production cells plus
stationary controls added on 31 October 2026. The final stationary grid has three
theta values, five sizes, 64 interface runs per input/cell and 15/15 passing binary
covariance checks at the unchanged 0.005 tolerance. Original failures are retained.

## Analyze the versioned repository data

From the repository root, with Python and requirements.txt installed:

    python -m pip install -r requirements.txt pytest
    python -m pytest -q
    python code/innovation_analysis.py --source data/raw/innovation_20261031_v3 --binary-source data/raw/innovation_binary_refined_20261031 --sizes 32 64 128 256 512 --out data/analysis/innovation_20261031_reproduced

The three innovation raw directories contain 15 NPZ archives each. The final
frozen analysis is results_refined/summary.json in this documentation directory.
The original 105-cell production protocol remains in
docs/submission_20261030/REPRODUCE.md. These are separate experimental grids.

## Independently verify the numerical supplement

Unpack manuscript/innovation_controls_data_20261031.zip into a separate empty
directory, enter that directory and run:

    python code/independent_supplement_audit.py
    python code/verify_submission_data.py

The independent audit needs NumPy and the Python standard library only. It imports
no project generator or analysis functions. It checks all payload hashes, recomputes
the original failures and final covariance gates, independently fits slopes and
bootstrap intervals, checks four-sign moments, reconstructs fitted spectra and
regenerates the first calibration record from its explicit seed for each refined
cell. It writes INDEPENDENT_AUDIT.json. The second command verifies the analysis
through the supplied project functions. DATA_INDEX.json and SEEDS.json identify
all relevant archive fields and per-run seeds.

Full reanalysis and paper Figure 7 reproduction from the same extraction:

    python code/innovation_analysis.py --source data/raw/innovation_20261031_v3 --binary-source data/raw/innovation_binary_refined_20261031 --sizes 32 64 128 256 512 --out reproduced_analysis
    python code/update_innovation_manuscript.py --summary summary.json --figure-only

Both independent and supplied analyses have been run from an actual extraction;
all fields of summary.json reproduce. See Supplemental_Material_核验报告.md and
supplement_actual_zip_audit.json for the archive hash and scope of those checks.

## Regenerate the refinement and validation

In a separate supplement extraction, after verifying the archived files:

    python code/refine_binary_controls.py
    python code/validate_large_binary_control.py

The first command uses 128 fresh 64-site calibration records and 64 fresh 64-site
validation records for theta=0.40 at L=32 and 64. Its seed root is 20261103, with
separate domains for training, Gaussian and uniform validation. The original
simultaneous bounds 0.010739/0.005850 become 0.003706/0.002553.

The second keeps the L=512 fitted spectrum and original interface records fixed
and adds a separate 128-run 64-site covariance-only cohort, seed root 20261105.
Its bound changes from 0.005310 to 0.002503. The tolerance remains 0.005.
These commands regenerate files and metadata; keep the untouched archive for
byte-level checks and compare regenerated numerical arrays with appropriate
floating-point tolerances. Original failed records remain under
innovation_binary_20261031_v3.

## Regenerate source campaigns

Full simulations are slower than reanalysis. Use fresh output directories rather
than overwriting archived evidence. Exact cache checks include producer source
hashes; a changed source file intentionally invalidates a cached result. Frozen
source snapshots retain the exact earlier numerical producers.

    python code/innovation_controls.py --thetas 0.1 0.25 0.4 --sizes 32 64 128 256 512 --runs 64 --workers 4 --out data/raw/reproduced_source

Baseline surrogate generation accepts an absolute --source directory and a fresh
--out directory through code/binary_covariance_control.py. The supplement contains
all final records, seeds and scripts needed to inspect the refinement without
rerunning the full campaigns.

## Manuscript and scope

The current PDF is manuscript/main.pdf (also BD_noise_controls_updated.pdf).
The editable package manuscript/overleaf_bd_paper.zip can be compiled on Overleaf
with the usual REVTeX/pdfLaTeX toolchain. The combined
PRE_submission_materials_20261031.zip also includes the numerical supplement and
README.TXT. No journal submission is performed by these reproduction commands.

Bootstrap units are independent runs, not stored times. Six-variant coupling is
retained within size; sizes are resampled independently. Surrogate comparisons
use held-out uniform records and are conditional on fitted spectra. Finite-lag
matching, finite-size slopes and non-detections are not claims of exact process
equivalence, asymptotic exponents or a universality-class change.
