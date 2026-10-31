import numpy as np
from bd_core import _u01
from innovation_controls import matched_weights, periodic_ar
import innovation_streaming as streaming


def test_streaming_pair_matches_explicit_circular_filter_and_independent_keys():
    w = np.array([.8, .3]); r = np.array([.2, .8])
    p, sites, seed = 128, 5, 123
    got = streaming.uniform_pair(w, r, p, sites, seed)
    sums = np.zeros((2, sites, p))
    v = matched_weights(w, r)
    base = (seed*0x9E3779B97F4A7C15) % 2**64
    for site in range(sites):
        for n in range(len(w)):
            z = np.array([_u01(np.uint64((base+(site*len(w)+n)*p+t) % 2**64))-.5 for t in range(p)])[None]
            z -= z.mean()
            x = periodic_ar(z, r[n])[0]
            sums[0, site] += w[n]*x
            sums[1, site] += v[n]*x
    np.testing.assert_array_equal(got, (sums > 0).transpose(0, 2, 1))


def test_formal_uniform_run_counter_intervals_are_disjoint():
    from bd_core import ffgn_n_components, smooth_ceil
    for theta in [.1, .25, .4]:
        for size in [32, 64, 128, 256, 512]:
            period = 2*smooth_ceil(round(15*size**1.5))
            length = size*period*ffgn_n_components(period)
            bases = []
            for run in range(64):
                children = np.random.SeedSequence([20261031, int(round(theta*1e6)), size, run]).spawn(2)
                seed = int(children[1].generate_state(1)[0])
                bases.append((seed*0x9E3779B97F4A7C15) % 2**64)
            ordered = sorted(bases)
            gaps = np.diff(np.array(ordered+[ordered[0]+2**64], dtype=object))
            assert min(gaps) > length
