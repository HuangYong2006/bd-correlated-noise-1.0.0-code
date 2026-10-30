"""Regression checks for independent-run diagnostics and honest uncertainty."""

import csv
import json
import sys

import numpy as np
import pytest

from diagnostics import (IForestDetector, discrete_hellinger_excess,
                         iforest_excess, onset_with_error)


def test_hellinger_error_is_outer_run_uncertainty_not_split_scatter():
    # Eight independent one-hot runs. Enumerating the Binomial(8, 1/2)
    # outer-bootstrap count and Hypergeometric(8, m, 4) inner halves gives
    # sd(E[H | outer ensemble]) = 0.07676388065406728. The old split sd is .122.
    ref = np.tile([0.5, 0.5], (8, 1))
    theta = np.array([[1.0, 0.0]] * 4 + [[0.0, 1.0]] * 4)
    result = discrete_hellinger_excess(ref, theta, n_splits=120, seed=0, n_boot=1000)
    assert result["dH_err"] == pytest.approx(0.07676388065406728, rel=0.20)


def test_zero_null_scatter_does_not_create_infinite_hellinger_z():
    ref = np.tile([0.5, 0.5], (8, 1))
    theta = np.tile([0.9, 0.1], (8, 1))
    result = discrete_hellinger_excess(ref, theta, n_splits=12, seed=0)
    assert np.isnan(result["z"])


def test_forest_has_disjoint_training_calibration_and_null_runs():
    ref = np.random.default_rng(17).normal(size=(60, 4))
    detector = IForestDetector(n_estimators=12, seed=3).fit(ref)
    assert getattr(detector, "null_scores", np.array([])).size == 12
    groups = detector.reference_indices
    assert len(set(groups["train"]) & set(groups["calibration"])) == 0
    assert len(set(groups["train"]) & set(groups["null"])) == 0
    assert len(set(groups["calibration"]) & set(groups["null"])) == 0
    # A finite held-out sample never establishes an exactly zero FPR.
    assert detector.phi_null_ci[1] > 0


def test_constant_forest_scores_have_undefined_standardized_effect():
    ref = np.ones((40, 4))
    detector = IForestDetector(n_estimators=8, seed=1).fit(ref)
    result = detector.score(ref[:8])
    assert np.isnan(result["z_mean"])


def test_repeated_snapshots_do_not_multiply_independent_forest_runs():
    rng = np.random.default_rng(8)
    ref = rng.normal(size=(40, 4))
    theta = rng.normal(loc=0.3, size=(40, 4))
    expected = iforest_excess(ref, theta, n_estimators=8, n_boot=3, seed=1)
    repeated_ref = np.repeat(ref[:, None, :], 4, axis=1)
    repeated_theta = np.repeat(theta[:, None, :], 4, axis=1)
    actual = iforest_excess(repeated_ref, repeated_theta,
                            n_estimators=8, n_boot=3, seed=1)
    assert actual["n_ref_runs"] == 40
    assert actual["dPhi"] == expected["dPhi"]
    assert actual["z_anom"] == expected["z_anom"]


def test_too_few_reference_runs_return_unavailable_forest_statistics():
    ref = np.arange(16, dtype=float).reshape(4, 4)
    result = iforest_excess(ref, ref, n_estimators=8, n_boot=2, seed=1)
    assert result["status"] == "insufficient_reference_runs"
    assert np.isnan(result["z_anom"])
    assert np.isnan(result["dPhi_err"])


def test_reference_forest_baseline_does_not_score_its_training_runs():
    ref = np.random.default_rng(24).normal(size=(60, 4))
    result = iforest_excess(ref, ref, n_estimators=12, n_boot=3, seed=2)
    assert result["dPhi"] == 0.0
    assert result["z_anom"] == 0.0


def test_missing_onset_errors_do_not_turn_into_zero_uncertainty():
    result = onset_with_error([0, 0.1, 0.2], [0, 1, 2],
                              [np.nan] * 3, 0.5)
    assert result["theta_star"] == pytest.approx(0.05)
    assert np.isnan(result["theta_star_err"])


def test_analysis_and_ablation_keep_snapshot_count_out_of_run_count(tmp_path, monkeypatch, capsys):
    import analyze
    import feature_ablation

    raw = tmp_path / "data" / "raw" / "unit"
    raw.mkdir(parents=True)
    rng = np.random.default_rng(13)
    for theta, token in [(0.0, "p0d0000"), (0.25, "p0d2500")]:
        # Binary heights guarantee the production model's {-1,0,1} increment
        # support while supplying four distinct late snapshots for twelve runs.
        profiles = rng.integers(0, 2, (12, 4, 32), dtype=np.int32)
        np.savez(raw / f"bd_{token}_L32.npz",
                 spec=json.dumps(dict(theta=theta, L=32, T=64, n_runs=12)),
                 N_ffgn=0, wall_time=0., check_times=np.array([16, 32, 48, 64]),
                 profiles=profiles, w_times=np.arange(1, 65),
                 w_series=rng.uniform(.2, .4, (12, 64)))
    monkeypatch.setattr(analyze, "ROOT", str(tmp_path))
    monkeypatch.setattr(sys, "argv", ["analyze.py", "--tag", "unit", "--splits", "4",
                                     "--hell-boot", "4", "--if-boot", "2", "--if-trees", "4"])
    analyze.main()
    with open(tmp_path / "data/analysis/unit_bd/diagnostics.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert {int(row["if_n_ref_runs"]) for row in rows} == {12}
    assert {int(row["if_n_null"]) for row in rows} == {3}
    assert {int(row["hell_n_bootstrap"]) for row in rows} == {4}

    monkeypatch.setattr(feature_ablation, "ROOT", str(tmp_path))
    monkeypatch.setattr(sys, "argv", ["feature_ablation.py", "--tag", "unit",
                                     "--if-boot", "2", "--if-trees", "4"])
    feature_ablation.main()
    result = json.loads((tmp_path / "data/analysis/ablation_unit.json").read_text())
    assert {cell["n_rows"] for cell in result["cells"]} == {12}
    assert result["statistics"]["unit"] == "run"
