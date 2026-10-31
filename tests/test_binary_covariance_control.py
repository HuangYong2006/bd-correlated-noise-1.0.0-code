import numpy as np
import pytest

import binary_covariance_control as binary
from innovation_controls import gaussian_sign_covariance


def test_positive_spectrum_fit_recovers_known_binary_covariance():
    period = 2048
    lag = np.unique(np.r_[np.arange(65), 128, 256, 512])
    spectra = binary.spectral_basis(period)
    target_s = .6*spectra[12]+.4*spectra[30]
    target_c = gaussian_sign_covariance(target_s, period)[lag]
    fitted, metadata = binary.fit_binary_spectrum(target_c, lag, period)
    assert (fitted >= 0).all()
    assert fitted[0] == 0
    assert metadata['max_training_residual'] < 1e-5
    np.testing.assert_allclose(gaussian_sign_covariance(fitted, period)[lag],
                               target_c, atol=1e-5)


def test_binary_match_rejects_invalid_target():
    with pytest.raises(ValueError):
        binary.fit_binary_spectrum(np.array([1., 1.1]), np.array([0, 1]), 128)
