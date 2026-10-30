"""Checks for the revised manuscript's width estimator and provenance guard."""
import numpy as np
import pytest
import revision_analysis as analysis


def subject(name):
    function = getattr(analysis, name, None)
    assert callable(function), f"revision analysis must implement {name}"
    return function


def test_late_width_samples_keep_runs_and_average_squared_width():
    w2 = np.array([np.arange(1., 11.), np.arange(11., 21.)])
    tail, previous = subject('late_width_samples')(w2)
    np.testing.assert_allclose(tail, [9.5, 19.5])
    np.testing.assert_allclose(previous, [7.5, 17.5])


def test_size_exponent_recovers_known_scaling_without_fit_residual_errors():
    squared = {L: np.full(8, L ** 1.5) for L in [16, 32, 64, 128]}
    result = subject('fit_size_exponent')(squared, n_boot=40, seed=7)
    assert result['alpha'] == pytest.approx(.75)
    assert result['alpha_err'] < 1e-12
    assert len(result['alpha_pairs']) == 3
    assert all(v['alpha'] == pytest.approx(.75) for v in result['alpha_pairs'])


def test_spectral_analysis_rejects_unversioned_old_data():
    check = subject('require_corrected_layout')
    with pytest.raises(ValueError, match='noise_layout'):
        check({'noise_layout': 'unversioned'}, 'fgn')
    check({'noise_layout': 'spectral-run-time-site-v2'}, 'fgn')
    check({'noise_layout': 'unversioned'}, 'shared')
