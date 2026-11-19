"""Paired spectral perturbations and joint convergence for theta=.25."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy import fft
from numba import set_num_threads

from bd_core import ffgn_n_components, smooth_ceil
from control_archives import atomic_savez
from controlled_roughening import growth_observables, noise_statistics, run_cell
from innovation_controls import gaussian_sign_covariance

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/raw/strict_binary_20261117'


def residual(sites):
    target = OUT/'residual'/f'L{sites}.npz'
    if target.exists():
        raise FileExistsError(target)
    source = ROOT/f'data/raw/innovation_binary_refined_20261031/theta0.25_L{sites}.npz'
    historical = ROOT/f'data/raw/innovation_20261031_v3/theta0.25_L{sites}.npz'
    with np.load(historical, allow_pickle=False) as z:
        spec = json.loads(str(z['spec']))
    with np.load(source, allow_pickle=False) as z:
        baseline = z['fitted_spectrum'].copy()
        lags = z['validation_lags'].copy()
    period = 2*spec['T']
    freq = fft.rfftfreq(period)
    directions = [np.ones(len(freq)), 1/(1+(period*freq/4)**2)]
    theory = gaussian_sign_covariance(baseline, period)
    spectra, mixtures = [baseline], []
    for direction in directions:
        direction[0] = 0
        direction *= fft.irfft(baseline, n=period)[0]/fft.irfft(direction, n=period)[0]
        lo, hi = 0., 1.
        for _ in range(35):
            f = (lo+hi)/2
            trial = (1-f)*baseline+f*direction
            dist = np.max(abs(gaussian_sign_covariance(trial, period)-theory))
            if dist <= .005:
                lo = f
            else:
                hi = f
        spectra.append((1-lo)*baseline+lo*direction)
        mixtures.append(lo)
    c = np.array([gaussian_sign_covariance(s, period) for s in spectra])
    np.testing.assert_allclose(np.max(abs(c[1:]-c[0]), axis=1), .005, atol=1e-9)
    amplitude = np.sqrt(spectra)
    def run(i):
        set_num_threads(1)
        rng = np.random.default_rng(np.random.SeedSequence([2026111711,sites,i]))
        eta = np.empty((3,period,sites), dtype=np.uint8)
        noise = None
        for first in range(0,sites,8):
            white = fft.rfft(rng.standard_normal((8,period)), axis=-1, workers=1)
            rows = fft.irfft(white[None]*amplitude[:,None,:], n=period, axis=-1, workers=1)
            eta[:,:,first:first+8] = (rows>0).transpose(0,2,1)
            if first == 0:
                noise = [noise_statistics(row,lags,[1]) for row in rows]
        w = [growth_observables(f,[spec['T'],period],spec['T']//2000)['late_w2'] for f in eta]
        return dict(late_w2=w, covariance=[n['sign_covariance'] for n in noise],
                    four_moment=np.concatenate([n['four_moment'] for n in noise]),
                    four_cumulant=np.concatenate([n['four_cumulant'] for n in noise]))
    with ThreadPoolExecutor(max_workers=2) as pool:
        records = []
        for i, row in enumerate(pool.map(run, range(64)), 1):
            records.append(row)
            if i%16 == 0:
                print(f'Residual L={sites}: {i}/64', flush=True)
    spec_out = dict(L=sites, T=spec['T'], period=period, theta=.25, seed=2026111711,
                    runs=64, variants=['baseline','white_mixture','low_frequency_mixture'],
                    diagnostic_sites=8, mixture_fractions=mixtures,
                    reference_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    atomic_savez(target, spec=json.dumps(spec_out), lags=lags, spectra=np.array(spectra),
                 theoretical_sign=c, **{k:np.stack([row[k] for row in records]) for k in records[0]})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('experiment', choices=['residual','joint'])
    args = parser.parse_args()
    if args.experiment == 'residual':
        (OUT/'residual').mkdir(parents=True, exist_ok=True)
        for sites in [64,128,256]:
            residual(sites)
    else:
        for j in range(3):
            for sites in [64,128]:
                base = smooth_ceil(round(15*sites**1.5))
                run_cell(.25, sites, 32, base, 2026111715+j, 2, OUT/f'joint{j}',
                         period=base*4*2**j, observation_times=[base*2**k for k in range(j+2)],
                         n_modes=ffgn_n_components(2*base)+4*j)


if __name__ == '__main__':
    main()

