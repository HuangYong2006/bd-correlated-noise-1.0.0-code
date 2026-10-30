"""Deterministic covariance checks, independent of a Monte Carlo noise floor."""

import numpy as np
import pytest

import check_generator_amplitude as amplitude


@pytest.mark.parametrize("shared", [False, True])
def test_single_ar_mode_matches_closed_form(shared):
    # X_t = .5 X_{t-1} + z_t, Var(z)=3, eta=2 X gives C(tau)=16*.5**tau.
    got = amplitude.ar_sum_covariance([2.0], [0.5], [0, 1, 2],
                                      shared=shared, innovation_variance=3.0)
    np.testing.assert_allclose(got, [16.0, 8.0, 4.0], rtol=1e-14)


@pytest.mark.parametrize("shared, expected", [
    (False, [28 / 5, 26 / 15, 3 / 5]),
    (True, [356 / 35, 362 / 105, 46 / 35]),
])
def test_two_ar_modes_keep_cross_terms_only_for_shared_innovations(shared, expected):
    # Diagonal terms give [28/5,26/15,3/5]; cross terms add [32/7,12/7,5/7].
    got = amplitude.ar_sum_covariance([1.0, 2.0], [0.5, 0.25], [0, 1, 2],
                                      shared=shared, innovation_variance=1.0)
    np.testing.assert_allclose(got, expected, rtol=1e-14)


def test_shared_covariance_matches_truncated_impulse_response_sum():
    weights = np.array([0.8, -0.3, 0.2])
    rates = np.array([0.2, 0.4, 0.6])
    impulse = np.array([np.dot(weights, rates ** j) for j in range(100)])
    lags = np.array([0, 1, 4, 9])
    expected = np.array([np.dot(impulse[:100 - lag], impulse[lag:])
                         for lag in lags]) / 12.0
    got = amplitude.ar_sum_covariance(weights, rates, lags)
    np.testing.assert_allclose(got, expected, rtol=1e-13, atol=1e-16)


@pytest.mark.parametrize("shared", [False, True])
def test_normalised_covariance_is_invariant_to_common_weight_scale(shared):
    weights = np.array([0.8, 0.3, 0.2])
    rates = np.array([0.2, 0.4, 0.6])
    lags = np.array([0, 1, 4, 9])
    first = amplitude.ar_sum_covariance(weights, rates, lags, shared=shared)
    scaled = amplitude.ar_sum_covariance(-7.0 * weights, rates, lags, shared=shared)
    np.testing.assert_allclose(scaled, 49.0 * first, rtol=1e-13)
    np.testing.assert_allclose(scaled / scaled[0], first / first[0], rtol=1e-13)


def test_amplitude_fit_does_not_rectify_a_signed_noise_floor():
    ac = np.array([1.0, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01])
    assert np.isnan(amplitude.fit_amplitude(ac, theta=0.1, lo=1, hi=6))


def test_amplitude_fit_preserves_a_positive_power_law():
    tau = np.arange(1, 21)
    ac = np.r_[1.0, 0.7 * tau ** -0.6]
    assert amplitude.fit_amplitude(ac, theta=0.2, lo=2, hi=20) == pytest.approx(0.7)
