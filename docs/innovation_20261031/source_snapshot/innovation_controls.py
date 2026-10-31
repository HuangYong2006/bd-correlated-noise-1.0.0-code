"""Stationary finite-record innovation controls; leaves production data intact.

All inputs use circular AR filters on a common length-2T grid with DC removed.
Gaussian mode sums can be sampled by a scalar spectral factor exactly in law.
Uniform independent sums explicitly generate one innovation stream per mode.
Raw covariance matching does NOT imply binary matching for uniform innovations.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy import fft, signal
from numba import set_num_threads

from bd_core import ffgn_n_components, ffgn_weights, smooth_ceil, _simulate_from_field
from innovation_streaming import uniform_pair

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ('gaussian_shared', 'gaussian_independent', 'gaussian_matched',
            'uniform_shared', 'uniform_independent', 'uniform_matched')
VERSION = 'innovation-circular-dc0-v3-disjoint-splitmix'


def matched_weights(weights, rates):
    """Positive AR weights matching the shared infinite stationary covariance."""
    a = weights * (weights[None, :] / (1-rates[:, None]*rates[None, :])).sum(axis=1)
    return np.sqrt(a * (1-rates**2))


def ar_spectra(weights, rates, period):
    phase = np.exp(-2j*np.pi*np.fft.rfftfreq(period))
    h = 1/(1-rates[:, None]*phase)
    shared = abs(weights @ h)**2
    independent = weights**2 @ abs(h)**2
    matched = matched_weights(weights, rates)**2 @ abs(h)**2
    for s in [shared, independent, matched]:
        s[0] = 0
    return shared, independent, matched


def gaussian_sign_covariance(spectrum, period):
    c = fft.irfft(spectrum, n=period)
    return 2/np.pi*np.arcsin(np.clip(c/c[0], -1, 1))


def periodic_ar(z, rate):
    """Solve x[t]=rate*x[t-1]+z[t] with x[-1]=x[P-1]."""
    out = signal.lfilter([1.], [1., -rate], z, axis=-1)
    p = z.shape[-1]
    correction = rate**np.arange(1, p+1) / (1-rate**p)
    return out + out[:, -1, None]*correction


@lru_cache(maxsize=1)
def spectral_filters(weight_tuple, rate_tuple, period):
    weights, rates = np.array(weight_tuple), np.array(rate_tuple)
    phase = np.exp(-2j*np.pi*np.fft.rfftfreq(period))
    transfer = sum(w/(1-r*phase) for w, r in zip(weights, rates))
    _, independent, matched = ar_spectra(weights, rates, period)
    transfer[0] = 0
    return transfer, np.sqrt(independent), np.sqrt(matched)


def independent_uniform_pair(weights, rates, period, sites, seed):
    """Common mode innovations; return original and covariance-matched weights."""
    return uniform_pair(weights, rates, period, sites, seed)


def noise_field(weights, rates, period, sites, seed, variant):
    """Binary (time, site) array; each site has an independent temporal record."""
    if variant not in VARIANTS:
        raise ValueError(variant)
    if variant in ['uniform_independent', 'uniform_matched']:
        pair = independent_uniform_pair(weights, rates, period, sites, seed)
        return pair[0 if variant == 'uniform_independent' else 1]
    rng = np.random.default_rng(seed)
    transfer, independent, matched = spectral_filters(tuple(weights), tuple(rates), period)
    if variant.startswith('gaussian'):
        amp = {'gaussian_shared': transfer,
               'gaussian_independent': independent,
               'gaussian_matched': matched}[variant]
    else:
        amp = transfer
    out = np.empty((period, sites), dtype=np.uint8)
    for start in range(0, sites, 8):
        n = min(8, sites-start)
        z = (rng.standard_normal((n, period))/np.sqrt(12)
             if variant.startswith('gaussian') else rng.uniform(-.5, .5, (n, period)))
        x = fft.irfft(fft.rfft(z, axis=-1)*amp, n=period, axis=-1)
        out[:, start:start+n] = (x > 0).T
    return out


def grow(eta, base_t):
    """Use one record for prefix-duration comparisons; widths are E_x[(h-hbar)^2]."""
    period, sites = eta.shape
    if period != 2*base_t:
        raise ValueError('record length must be twice base duration')
    stride = max(1, base_t//2000)
    times = np.arange(stride, period+1, stride)
    profiles = np.zeros((1, 2, sites), dtype=np.int32)
    widths = np.zeros((1, len(times)), dtype=np.float32)
    _simulate_from_field(eta[None], 0, sites, period, stride,
                         np.array([base_t, period]), profiles, widths)
    if np.max(abs(profiles-np.roll(profiles, 1, axis=-1))) > 1:
        raise AssertionError('parallel BD single-step invariant failed')
    late = [widths[0, (times > int(.8*t)) & (times <= t)].mean(dtype=float)
            for t in [base_t, period]]
    centered = profiles[0].astype(float)-profiles[0].mean(axis=-1, keepdims=True)
    # Full FFT: sum_k S(k) equals spatial squared width (Parseval).
    structure = abs(fft.fft(centered, axis=-1))**2/sites**2
    separations = np.unique([1, 2, 4, max(1, sites//8), max(1, sites//4)])
    height_difference = np.array([(centered-np.roll(centered, int(d), axis=-1))**2
                                 for d in separations]).mean(axis=-1).T
    return dict(late_w2=np.array(late), profiles=profiles[0],
                structure_factor=structure, separations=separations,
                height_difference=height_difference,
                flat_fraction=(np.diff(np.c_[profiles[0], profiles[0, :, :1]], axis=-1)==0).mean(axis=-1))


def noise_diagnostics(eta, lags):
    """Raw sign moments without record-wise demeaning; independent sites pooled."""
    x = 2*eta[:, :min(16, eta.shape[1])].T.astype(float)-1
    power = (abs(fft.rfft(x, axis=-1))**2).mean(axis=0)
    ac = fft.irfft(power, n=x.shape[1])/x.shape[1]
    # Persistence is P(all positive for ell successive times), not a fitted tail.
    lengths = np.array([1, 2, 4, 8, 16, 32])
    persistence = []
    for length in lengths:
        sums = np.cumsum(np.c_[np.zeros(len(x)), x > 0], axis=-1)
        persistence.append(((sums[:, length:]-sums[:, :-length]) == length).mean())
    return dict(sign_covariance=ac[lags], deposition_mean=float(eta.mean()),
                spatial_sign_correlation=float((x*np.roll(x, 1, axis=0)).mean()),
                persistence=np.array(persistence))


def one_run(spec, run):
    set_num_threads(1)
    weights, rates = ffgn_weights(spec['theta'], spec['N'])
    seed = np.random.SeedSequence([spec['seed'], int(round(spec['theta']*1e6)), spec['L'], run])
    # Common base seeds couple Gaussian variants and independently weighted
    # uniform variants. Whole-run pairing is preserved in inference.
    gaussian_seed, uniform_seed = [int(s.generate_state(1)[0]) for s in seed.spawn(2)]
    rows = []
    for variant in VARIANTS[:4]:
        eta = noise_field(weights, rates, 2*spec['T'], spec['L'],
                          gaussian_seed if variant.startswith('gaussian') else uniform_seed,
                          variant)
        rows.append({**grow(eta, spec['T']), **noise_diagnostics(eta, spec['lags'])})
    del eta
    pair = independent_uniform_pair(weights, rates, 2*spec['T'], spec['L'], uniform_seed)
    for eta in pair:
        rows.append({**grow(eta, spec['T']), **noise_diagnostics(eta, spec['lags'])})
    return rows


def run_cell(theta, sites, runs, base_t, seed, workers, outdir):
    spec = dict(version=VERSION, theta=theta, L=sites, T=base_t, runs=runs,
                N=ffgn_n_components(2*base_t), seed=seed,
                lags=np.unique(np.r_[0, 1, 2, 4, 8, 16, 32, 64,
                                    np.geomspace(128, max(128, base_t//2), 12).astype(int)]).tolist(),
                variants=list(VARIANTS), record_length=2*base_t,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                streaming_sha256=hashlib.sha256((ROOT/'code/innovation_streaming.py').read_bytes()).hexdigest(),
                solver_sha256=hashlib.sha256((ROOT/'code/bd_core.py').read_bytes()).hexdigest())
    path = outdir/f'theta{theta:.2f}_L{sites}.npz'
    if path.exists():
        with np.load(path, allow_pickle=False) as z:
            if json.loads(str(z['spec'])) != spec:
                raise ValueError(f'incompatible existing cell: {path}')
        print(f'Validated cached {path.name}', flush=True)
        return path
    t0 = time.time()
    weights, rates = ffgn_weights(theta, spec['N'])
    spectral_filters(tuple(weights), tuple(rates), 2*base_t)
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, result in enumerate(pool.map(lambda run: one_run(spec, run), range(runs))):
            rows.append(result)
            if (i+1) % 8 == 0 or i+1 == runs:
                print(f'theta={theta:.2f} L={sites} {i+1}/{runs} runs, {time.time()-t0:.1f}s', flush=True)
    data = {key: np.array([[row[v][key] for v in range(len(VARIANTS))] for row in rows])
            for key in rows[0][0]}
    weights, rates = ffgn_weights(theta, spec['N'])
    spectra = ar_spectra(weights, rates, 2*base_t)
    targets = np.array([gaussian_sign_covariance(s, 2*base_t)[spec['lags']] for s in spectra])
    np.savez_compressed(path, spec=json.dumps(spec), **data,
                        gaussian_sign_targets=targets, wall_seconds=time.time()-t0)
    print(f'Saved {path.name}: {time.time()-t0:.1f}s', flush=True)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--thetas', nargs='+', type=float, default=[.1, .25, .4])
    ap.add_argument('--sizes', nargs='+', type=int, default=[32, 64, 128, 256, 512])
    ap.add_argument('--runs', type=int, default=64)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--seed', type=int, default=20261031)
    ap.add_argument('--pilot', action='store_true', help='use T=512; never mix with formal cells')
    ap.add_argument('--out', type=Path, default=ROOT/'data/raw/innovation_20261031_v3')
    args = ap.parse_args()
    if args.runs < 2 or args.workers < 1 or any(s < 8 for s in args.sizes):
        ap.error('need >=2 runs, >=1 worker and L>=8')
    if any(not 0 <= theta < .5 for theta in args.thetas):
        ap.error('theta must be in [0, .5)')
    if args.pilot:
        args.out = args.out/'pilot'
    args.out.mkdir(parents=True, exist_ok=True)
    for sites in args.sizes:
        for theta in args.thetas:
            base_t = 512 if args.pilot else smooth_ceil(round(15*sites**1.5))
            run_cell(theta, sites, args.runs, base_t, args.seed, args.workers, args.out)


if __name__ == '__main__':
    main()
