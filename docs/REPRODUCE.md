# Reproducing the numerical data

Start with [the current numerical reproduction instructions](strict_binary_20261117/REPRODUCE.md). Download the versioned numerical archive from the repository Release and extract it into an empty directory. It contains the raw arrays and frozen inputs that are deliberately not stored directly in Git.

The current instructions specify parity, covariance perturbation and joint-limit experiments, seeds, resampling and source-data files. The [historical controlled-experiment instructions](research_refocus_20261117/REPRODUCE.md) give earlier observation-time, period and mode controls.

Use the supplied requirements and the documented environment. From the extracted root:

```powershell
python code/strict_binary_analysis.py
python code/refocus_figures.py --out figs/strict_binary_20261117 --notes docs/strict_binary_20261117
python code/strict_binary_figure.py
python -m pytest -q -p no:cacheprovider
```

For a complete historical calibration audit, extract only missing files from the nested innovation archive and retain newer code/manifests. The package README explains coverage limits. This repository and Release distribute numerical data and reproducibility materials, not manuscript drafts or article-source bundles.
