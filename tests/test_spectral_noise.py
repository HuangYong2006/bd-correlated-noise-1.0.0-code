"""Production spectral fields must correlate in time, independently per site."""

import json

import numpy as np
import pytest

from bd_core import (SimSpec, generate_spectral_eta, load_result, save_result,
                     simulate_spectral)


def _field(spectrum, L, T, n_runs, seed=20261026, rows_per_chunk=5,
           runs_per_chunk=2):
    chunks = list(generate_spectral_eta(
        0.25, L, T, n_runs, seed, spectrum=spectrum,
        rows_per_chunk=rows_per_chunk, max_bytes=runs_per_chunk * T * L))
    expected_starts = list(range(0, n_runs, runs_per_chunk))
    assert [r0 for r0, _ in chunks] == expected_starts
    for r0, eta in chunks:
        assert eta.shape == (min(runs_per_chunk, n_runs-r0), T, L)
        assert eta.dtype == np.uint8
        assert np.all((eta == 0) | (eta == 1))
    return np.concatenate([eta for _, eta in chunks], axis=0)


def _correlation(a, b):
    a = a.astype(float) - a.mean()
    b = b.astype(float) - b.mean()
    return np.mean(a*b) / np.sqrt(np.mean(a*a)*np.mean(b*b))


@pytest.mark.parametrize("spectrum", ["powerlaw", "fgn"])
@pytest.mark.parametrize("L", [3, 17])
def test_production_field_correlates_in_time_not_space(spectrum, L):
    # Swapping the time/site layout leaves the 1-D filter correct but violates
    # both checks. Non-divisible row chunks cross site and run boundaries.
    eta = _field(spectrum, L, T=16384, n_runs=3)
    temporal = _correlation(eta[:, :-1, :], eta[:, 1:, :])
    spatial = _correlation(eta[:, :, :-1], eta[:, :, 1:])
    assert temporal > 0.15, f"temporal lag-1 correlation is {temporal}"
    assert abs(spatial) < 0.03, f"spatial lag-1 correlation is {spatial}"
    assert abs(_correlation(eta[0], eta[1])) < 0.03


@pytest.mark.parametrize("spectrum", ["powerlaw", "fgn"])
@pytest.mark.parametrize("L,n_runs", [(1, 1), (3, 5), (7, 4)])
def test_seeded_samples_are_invariant_to_row_and_run_chunks(spectrum, L, n_runs):
    expected = _field(spectrum, L, 256, n_runs, rows_per_chunk=1,
                      runs_per_chunk=1)
    for rows, runs in [(2, 2), (5, 3), (64, n_runs)]:
        actual = _field(spectrum, L, 256, n_runs, rows_per_chunk=rows,
                        runs_per_chunk=runs)
        np.testing.assert_array_equal(actual, expected)
    repeated = _field(spectrum, L, 256, n_runs, rows_per_chunk=1,
                      runs_per_chunk=1)
    np.testing.assert_array_equal(repeated, expected)
    different_seed = _field(spectrum, L, 256, n_runs, seed=43)
    assert not np.array_equal(different_seed, expected)


@pytest.mark.parametrize("spectrum", ["powerlaw", "fgn"])
def test_spectral_simulation_archive_keeps_layout_version(tmp_path, spectrum):
    spec = SimSpec(theta=0.25, L=5, T=64, n_runs=2, n_check=2,
                   n_w=16, seed=23, spectrum=spectrum)
    result = simulate_spectral(spec)
    assert result["noise_layout"] == "spectral-run-time-site-v2"
    archive = tmp_path / "corrected.npz"
    save_result(str(archive), result)
    loaded = load_result(str(archive))
    assert loaded["noise_layout"] == result["noise_layout"]
    assert isinstance(loaded["noise_layout"], str)
    np.testing.assert_array_equal(loaded["profiles"], result["profiles"])


def test_legacy_archive_is_not_labelled_as_corrected(tmp_path):
    archive = tmp_path / "legacy.npz"
    np.savez_compressed(archive, spec=json.dumps({"theta": 0.25}),
                        N_ffgn=0, wall_time=0.0,
                        profiles=np.zeros((1, 1, 3), dtype=np.int32))
    assert load_result(str(archive))["noise_layout"] == "unversioned"
