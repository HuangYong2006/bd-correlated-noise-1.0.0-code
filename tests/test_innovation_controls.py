"""Physical and numerical checks for the innovation-control campaign."""
import numpy as np
import pytest

import innovation_controls as controls
from check_generator_amplitude import ar_sum_covariance


def test_reweighted_independent_covariance_matches_shared_at_all_lags():
    w = np.array([0.8, 0.3, 0.2])
    r = np.array([0.2, 0.6, 0.95])
    v = controls.matched_weights(w, r)
    lags = np.arange(200)
    np.testing.assert_allclose(ar_sum_covariance(w, r, lags),
                               ar_sum_covariance(v, r, lags, shared=False), rtol=2e-14)
    shared, independent, matched = controls.ar_spectra(w, r, 1024)
    np.testing.assert_allclose(shared, matched, rtol=2e-14)
    assert shared[0] == independent[0] == matched[0] == 0
    assert not np.allclose(shared[1:], independent[1:])


def test_periodic_filter_solves_circular_ar_equation():
    z = np.random.default_rng(71).uniform(-0.5, 0.5, (3, 128))
    for r in [0., 0.8, 0.999999]:
        x = controls.periodic_ar(z, r)
        np.testing.assert_allclose(x - r*np.roll(x, 1, axis=-1), z, atol=1e-10)


def test_gaussian_matching_survives_thresholding_and_spatial_assembly():
    w = np.array([0.8, 0.3]); r = np.array([0.2, 0.8])
    shared, _, matched = controls.ar_spectra(w, r, 4096)
    expected = controls.gaussian_sign_covariance(shared, 4096)
    np.testing.assert_allclose(expected,
                               controls.gaussian_sign_covariance(matched, 4096), atol=1e-14)
    eta = controls.noise_field(w, r, 4096, 128, 123, 'gaussian_shared')
    assert eta.shape == (4096, 128)
    assert eta.dtype == np.uint8
    assert 0.49 < eta.mean() < 0.51
    sign = 2*eta.astype(float)-1
    for lag in [1, 2, 4]:
        assert abs((sign[:-lag]*sign[lag:]).mean()-expected[lag]) < 0.012
    assert abs((sign*np.roll(sign, 1, axis=1)).mean()) < 0.012
    np.testing.assert_array_equal(eta, controls.noise_field(w, r, 4096, 128, 123,
                                                          'gaussian_shared'))


def test_growth_observables_match_direct_update():
    eta = np.random.default_rng(9).integers(0, 2, (96, 12), dtype=np.uint8)
    got = controls.grow(eta, 48)
    h = np.zeros(12, dtype=np.int64)
    widths = []
    for e in eta:
        h = np.maximum.reduce([h+e, np.roll(h, 1), np.roll(h, -1)])
        widths.append(h.var())
    np.testing.assert_array_equal(got['profiles'][-1], h)
    assert got['late_w2'][0] == pytest.approx(np.mean(widths[38:48]))
    assert got['late_w2'][1] == pytest.approx(np.mean(widths[76:96]))
    assert got['structure_factor'][-1].sum() == pytest.approx(h.var())


def test_uniform_variants_are_unbranched_at_theta_zero():
    from bd_core import ffgn_weights
    w, r = ffgn_weights(0., 16)
    for variant in ['uniform_shared', 'uniform_independent', 'uniform_matched']:
        eta = controls.noise_field(w, r, 2048, 32, 31, variant)
        sign = 2*eta.astype(float)-1
        assert (sign[:-1]*sign[1:]).mean() > 0.08
    with pytest.raises(ValueError):
        controls.noise_field(w, r, 10, 2, 1, 'typo')


def test_paired_uniform_generation_preserves_explicit_independent_streams():
    w = np.array([0.8, 0.3]); r = np.array([0.2, 0.8])
    pair = controls.independent_uniform_pair(w, r, 512, 12, 73)
    for i, name in enumerate(['uniform_independent', 'uniform_matched']):
        np.testing.assert_array_equal(pair[i], controls.noise_field(w, r, 512, 12, 73, name))
