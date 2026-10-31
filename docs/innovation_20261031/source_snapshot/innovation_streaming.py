"""Memory-bounded two-pass circular AR sums with independent counter streams."""
import numpy as np
from numba import njit
from bd_core import _u01


@njit(cache=True, nogil=True)
def uniform_pair(weights, rates, period, sites, seed):
    # Spread run seeds through the 64-bit counter space. Adding 32-bit seeds
    # directly would reuse long intervals across independent large runs.
    run_base = np.uint64(seed)*np.uint64(0x9E3779B97F4A7C15)
    modes = len(weights)
    matched = np.empty(modes)
    for n in range(modes):
        a = 0.
        for m in range(modes):
            a += weights[m]/(1-rates[n]*rates[m])
        matched[n] = np.sqrt(weights[n]*a*(1-rates[n]**2))
    state = np.zeros((sites, modes))
    means = np.zeros((sites, modes))
    # First pass computes the exact periodic initial state and innovation DC.
    for site in range(sites):
        for n in range(modes):
            total = 0.
            value = 0.
            base = run_base+np.uint64((site*modes+n)*period)
            for t in range(period):
                z = _u01(base+np.uint64(t))-.5
                total += z
                value = rates[n]*value+z
            means[site, n] = total/period
            state[site, n] = value/(1-rates[n]**period)-means[site, n]/(1-rates[n])
    out = np.empty((2, period, sites), dtype=np.uint8)
    for t in range(period):
        for site in range(sites):
            original = 0.
            corrected = 0.
            for n in range(modes):
                key = run_base+np.uint64((site*modes+n)*period+t)
                z = _u01(key)-.5-means[site, n]
                value = rates[n]*state[site, n]+z
                state[site, n] = value
                original += weights[n]*value
                corrected += matched[n]*value
            out[0, t, site] = original > 0
            out[1, t, site] = corrected > 0
    return out
