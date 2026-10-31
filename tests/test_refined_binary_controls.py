import numpy as np

from refine_binary_controls import all_site_diagnostics, fresh_seed


def test_refined_diagnostics_use_all_64_sites():
    eta = np.ones((128, 64), dtype=np.uint8)
    eta[::2, :16] = 0
    lags = np.array([0, 1, 2, 3])
    result = all_site_diagnostics(eta, lags)
    sign = 2*eta.astype(float)-1
    expected = [(sign*np.roll(sign, int(lag), axis=0)).mean() for lag in lags]
    np.testing.assert_allclose(result['sign_covariance'], expected, atol=1e-14)
    assert result['sign_covariance'][1] == .5


def test_refinement_training_and_validation_seed_domains_are_disjoint():
    seeds = [fresh_seed(size, domain, run) for size in [32, 64]
             for domain, count in [(0, 128), (1, 64), (2, 64)] for run in range(count)]
    assert len(set(seeds)) == len(seeds)
