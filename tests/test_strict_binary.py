from fractions import Fraction

import numpy as np

from strict_binary import binary_field, exact_width


def test_parity_field_pairwise_white_with_consecutive_four_sign_signal():
    rng = np.random.default_rng(1021)
    eta = binary_field(rng, 16384, 64, True)
    s = 2*eta.astype(np.int8)-1
    assert abs(s.mean()) < .005
    for lag in [1,2,3,4,19]:
        assert abs((s[:-lag]*s[lag:]).mean()) < .005
    assert abs((s[:,1:]*s[:,:-1]).mean()) < .005
    assert abs((s[:-3]*s[1:-2]*s[2:-1]*s[3:]).mean()-.25) < .005


def test_exact_width_effect_starts_at_fourth_sweep():
    result = exact_width(4)
    assert list(map(Fraction,result['difference'])) == [0,0,0,Fraction(2595,262144)]
