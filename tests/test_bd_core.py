"""
Unit tests for the core invariants of the BD model and its noise generators.

These are deliberately small and fast: they check *structure*, not physics. The
physics validation (alpha = 1/2, beta = 1/3 at theta = 0) is a campaign-scale
check and lives in code/validate_model.py.

Run from the repository root:

    pytest
    # or
    python -m pytest tests -q

The tests import the simulation modules the same way the drivers do, i.e. with
code/ on sys.path (configured via [tool.pytest.ini_options] pythonpath in
pyproject.toml).
"""

from __future__ import annotations

import numpy as np
import pytest

from bd_core import (SimSpec, ffgn_n_components, ffgn_weights, fgn_spectral_amp,
                     smooth_ceil, simulate)
from features import height_difference_correlation, snapshot_features
from diagnostics import (hellinger, increment_pmf, discrete_hellinger_excess,
                         onset_theta, terrace_length_pmf)


# --------------------------------------------------------------------------
# noise generators
# --------------------------------------------------------------------------

def test_ffgn_component_count_grows_with_T():
    """Mandelbrot's rule: more components are needed to resolve longer records."""
    assert ffgn_n_components(1 << 10) <= ffgn_n_components(1 << 20)
    # a floor of 16 components is imposed regardless of T
    assert ffgn_n_components(16) >= 16


def test_ffgn_amplitude_does_not_vanish_at_zero_theta():
    """The implemented weight norm remains finite as theta -> 0.

    This structural check does not measure a covariance-tail amplitude or
    establish a universality class. Shared AR innovations add cross terms.
    """
    W0, _ = ffgn_weights(0.0, 32)
    W1, _ = ffgn_weights(0.1, 32)
    a0 = float(np.sqrt(np.sum(W0 ** 2)))
    a1 = float(np.sqrt(np.sum(W1 ** 2)))
    assert a0 > 0.5, "FFGN amplitude at theta=0 should be O(1), not vanishing"
    # and it is nearly theta-independent
    assert abs(a1 - a0) / a0 < 0.5


def test_exact_fgn_amplitude_is_white_at_zero_theta():
    """The exact fGn spectral density must be flat at theta = 0.

    The full infinite alias sum is flat at theta=0; the implemented finite sum
    approximates it on nonzero modes. This is only a loose in-band check, not a
    claim that a finite record with DC removed is independent white noise.
    """
    T = 256
    amp = fgn_spectral_amp(0.0, T)
    band = amp[8:T // 2]
    assert band.size > 0
    assert band.min() > 0.0
    assert band.max() / band.min() < 5.0, "exact fGn spectrum is not flat at theta=0"


def test_smooth_ceil_is_five_smooth_and_ge():
    """smooth_ceil must return a 5-smooth integer >= n (fast FFT length)."""
    for n in (1, 7, 100, 1000, 21722):
        m = smooth_ceil(n)
        assert m >= n
        r = m
        for p in (2, 3, 5):
            while r % p == 0:
                r //= p
        assert r == 1, f"{m} is not 5-smooth (n={n})"


# --------------------------------------------------------------------------
# model structure
# --------------------------------------------------------------------------

def test_profiles_are_single_step_from_flat():
    """|Delta h| <= 1 must hold everywhere, from a flat initial condition.

    This is what reduces the increment pmf to a one-parameter family and makes
    'faceting' impossible to explain by overhangs.
    """
    spec = SimSpec(theta=0.0, L=64, T=200, n_runs=4, n_check=4, n_w=200, seed=7)
    res = simulate(spec)
    p = res["profiles"].astype(np.int64)
    d = np.abs(p - np.roll(p, -1, axis=-1))
    assert d.max() <= 1, "the BD rule violated |Delta h| <= 1"


def test_increment_pmf_sums_to_one_and_has_three_bins():
    spec = SimSpec(theta=0.0, L=64, T=100, n_runs=4, n_check=2, n_w=100, seed=11)
    res = simulate(spec)
    pmf = increment_pmf(res["profiles"][:, -1, :])
    assert pmf.shape == (4, 3)
    assert np.allclose(pmf.sum(axis=1), 1.0)


def test_up_and_down_rates_are_equal_for_bd():
    """sum_x Delta h = 0 telescopes exactly on a periodic lattice.

    With support on {-1, 0, +1} this forces P(+1) = P(-1) per configuration,
    hence skew(Delta h) = 0 to machine precision.
    """
    spec = SimSpec(theta=0.0, L=128, T=150, n_runs=4, n_check=2, n_w=150, seed=3)
    res = simulate(spec)
    p = res["profiles"][:, -1, :].astype(np.int64)
    delta = p - np.roll(p, -1, axis=-1)
    assert np.array_equal(delta.sum(axis=-1), np.zeros(p.shape[0], dtype=np.int64))
    pmf = increment_pmf(res["profiles"][:, -1, :])
    assert np.allclose(pmf[:, 0], pmf[:, 2], atol=1e-12)


def test_correlation_c1_is_positive_on_a_rough_interface():
    spec = SimSpec(theta=0.0, L=128, T=150, n_runs=4, n_check=2, n_w=150, seed=5)
    res = simulate(spec)
    C = height_difference_correlation(res["profiles"][:, -1, :].astype(np.float64))
    assert C.shape[-1] == 9
    # Only separations well below L/2 are meaningful: for r >= L/2 the periodic
    # estimator compares sites a full box apart and the correlation collapses to
    # zero by construction. R_LIST is (1,2,4,8,16,32,64,128,256), so at L=128 the
    # first seven entries (r <= 64) are the physical ones.
    assert np.all(C[:, :7] > 0)
    # C(r) must grow with r. A single short run is too noisy for a per-run
    # monotonicity claim, so we assert it on the ensemble average, which is the
    # quantity the analysis actually fits.
    Cmean = C[:, :7].mean(axis=0)
    assert np.all(np.diff(Cmean) > 0), f"C(r) not increasing: {Cmean}"


# --------------------------------------------------------------------------
# diagnostics
# --------------------------------------------------------------------------

def test_hellinger_bounds_and_identity():
    p = np.array([0.2, 0.3, 0.5])
    assert hellinger(p, p) == pytest.approx(0.0, abs=1e-12)
    # disjoint support saturates at 1
    assert hellinger(np.array([1.0, 0.0]), np.array([0.0, 1.0])) == pytest.approx(1.0)
    # symmetry
    q = np.array([0.1, 0.6, 0.3])
    assert hellinger(p, q) == pytest.approx(hellinger(q, p))


def test_terrace_pmf_sums_to_one():
    spec = SimSpec(theta=0.0, L=64, T=100, n_runs=3, n_check=2, n_w=100, seed=13)
    res = simulate(spec)
    pmf = terrace_length_pmf(res["profiles"][:, -1, :], max_len=16)
    assert pmf.shape == (3, 17)
    assert np.allclose(pmf.sum(axis=1), 1.0)


def test_split_half_null_is_centred_on_zero():
    """The split-half calibration must report no *excess* for identical data.

    Note the estimator is a difference of two finite-sample Hellinger distances,
    so it is centred on zero but not identically zero: each half-ensemble has its
    own sampling noise. The correct check is therefore that the excess is small
    compared with the scatter of the null floor itself, which is exactly the
    property that makes the calibration usable as a significance measure.
    """
    spec = SimSpec(theta=0.2, L=64, T=300, n_runs=32, n_check=1, n_w=300, seed=17)
    res = simulate(spec)
    pmf = increment_pmf(res["profiles"][:, -1, :])
    r = discrete_hellinger_excess(pmf, pmf, n_splits=40, seed=0)
    assert r["Hnull_mean"] > 0.0
    assert np.isfinite(r["z"])
    assert abs(r["dH"]) < 5.0 * max(r["Hnull_std"], 1e-12)
    assert abs(r["z"]) < 5.0


def test_onset_interpolation():
    thetas = np.array([0.0, 0.1, 0.2, 0.3])
    values = np.array([0.0, 0.0, 1.0, 1.0])
    # linear interpolation should land on the midpoint between 0.1 and 0.2
    assert onset_theta(thetas, values, 0.5) == pytest.approx(0.15)
    # no crossing -> nan, never a silent 0
    assert np.isnan(onset_theta(thetas, values, 5.0))


def test_onset_on_decreasing_series():
    thetas = np.array([0.0, 0.1, 0.2, 0.3])
    values = np.array([1.0, 1.0, 0.0, 0.0])
    assert onset_theta(thetas, values, 0.5, increasing=False) == pytest.approx(0.15)


# --------------------------------------------------------------------------
# features
# --------------------------------------------------------------------------

def test_feature_matrix_shape_and_finiteness():
    spec = SimSpec(theta=0.0, L=64, T=120, n_runs=4, n_check=2, n_w=120, seed=19)
    res = simulate(spec)
    names, F, extra = snapshot_features(res["profiles"])
    assert F.shape[0] == 4
    assert F.shape[-1] == len(names)
    # amplitude-free features must all be finite
    assert np.all(np.isfinite(F)), f"non-finite features: {names}"
    assert "alpha_s" in names
    assert "log_mean_flat_run" in names
    # the amplitude itself is carried in `extra`, deliberately outside the
    # feature vector, so that no feature can see the non-universal growth
    assert extra["W"].shape == (4, 2)
