"""Exact and sampled growth under strictly pairwise-matched binary inputs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
from numba import set_num_threads

from bd_core import smooth_ceil
from control_archives import atomic_savez
from controlled_roughening import growth_observables

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/raw/strict_binary_20261117'
LAGS = [1, 2, 3, 4, 8, 16, 32, 64, 128, 256]


def binary_field(rng, length, sites, parity):
    if not parity:
        return rng.integers(0, 2, (length, sites), dtype=np.uint8)
    # Independent stationary phases. Full blocks extend beyond both endpoints.
    phase = rng.integers(0, 4, sites)
    result = np.empty((length, sites), dtype=np.uint8)
    blocks = (length+6)//4
    for site, offset in enumerate(phase):
        bits = rng.integers(0, 2, (blocks, 4), dtype=np.uint8)
        bits[:, 3] = bits[:, 0] ^ bits[:, 1] ^ bits[:, 2]
        result[:, site] = bits.ravel()[offset:offset+length]
    return result


def exact_width(sites):
    """Enumerate all four-sweep inputs with exact integer probability weights."""
    first, second = [0]*4, [0]*4
    shifts = np.arange(4*sites, dtype=np.uint64)
    for start in range(0, 2**(4*sites), 65536):
        indices = np.arange(start, min(start+65536, 2**(4*sites)), dtype=np.uint64)
        bits = ((indices[:, None] >> shifts) & 1).astype(np.int64).reshape(-1, 4, sites)
        # Each site's four-bit marginal has weight (4 + product of signs)/64.
        numerator = np.prod(4+np.prod(2*bits-1, axis=1), axis=1)
        h = np.zeros((len(indices), sites), dtype=np.int64)
        for t in range(4):
            h = np.maximum.reduce((h+bits[:, t], np.roll(h, 1, axis=1), np.roll(h, -1, axis=1)))
            scaled = sites*(h*h).sum(axis=1)-h.sum(axis=1)**2
            first[t] += int(scaled.sum())
            second[t] += int((scaled*numerator).sum())
    a = [Fraction(x, sites**2*16**sites) for x in first]
    b = [Fraction(x, sites**2*64**sites) for x in second]
    assert a[:3] == b[:3]
    return dict(L=sites, iid=list(map(str,a)), parity=list(map(str,b)),
                difference=[str(y-x) for x,y in zip(a,b)],
                relative_width_change=(float(b[-1]/a[-1])**.5-1))


def moments(field):
    # Diagnostic sites in physical order, without circular wrap in time.
    s = 2*field[:, :16].astype(np.int8)-1
    corr = np.array([(s[:-lag]*s[lag:]).mean() for lag in LAGS])
    fourth = np.array([(s[:-3*d]*s[d:-2*d]*s[2*d:-d]*s[3*d:]).mean() for d in [1,2,4,8]])
    return dict(deposition=field.mean(), corr=corr, fourth=fourth,
                spatial=(s[:, :-1]*s[:, 1:]).mean())


def simulate(sites, runs, workers):
    base = smooth_ceil(round(15*sites**1.5))
    path = OUT/f'L{sites}.npz'
    if path.exists():
        raise FileExistsError(path)
    def run(i):
        set_num_threads(1)
        outputs = []
        for parity in [0, 1]:
            rng = np.random.default_rng(np.random.SeedSequence([2026111707, sites, i, parity]))
            field = binary_field(rng, 2*base, sites, parity)
            record = growth_observables(field, [base,2*base], max(1,base//1000))
            # Width snapshots at the first six times are sampled exactly.
            early = growth_observables(field[:16], [1,2,3,4,8,16])
            record['early_w2'] = early['profiles'].var(axis=1)
            outputs.append({**record, **moments(field)})
        return {key: np.stack([out[key] for out in outputs]) for key in outputs[0]}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        records = []
        for i, sample in enumerate(pool.map(run, range(runs)), 1):
            records.append(sample)
            if i%32 == 0:
                print(f'Parity L={sites}: {i}/{runs}', flush=True)
    data = {k: np.stack([r[k] for r in records]) for k in records[0]}
    assert np.all(np.abs(data['profiles']-np.roll(data['profiles'],1,axis=-1))<=1)
    np.testing.assert_allclose(data['structure_factor'].sum(axis=-1)/sites,
                               data['profiles'].var(axis=-1), atol=1e-9)
    metadata = dict(seed=2026111707, seed_fields=['seed','L','run','variant'], T=base, L=sites,
                    runs=runs, observation_times=[base,2*base], early_times=[1,2,3,4,8,16],
                    variants=['iid','parity'], covariance_lags=LAGS,
                    fourth_spacings=[1,2,4,8], diagnostic_sites=16,
                    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    solver_sha256=hashlib.sha256((ROOT/'code/bd_core.py').read_bytes()).hexdigest())
    atomic_savez(path, spec=json.dumps(metadata), **data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes', type=int, nargs='+', default=[64,128,256,512])
    parser.add_argument('--runs', type=int, default=128)
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    exact = [exact_width(L) for L in [3,4,5]]
    (OUT/'exact.json').write_text(json.dumps(exact, indent=2)+'\n')
    print('Exact:', exact, flush=True)
    for sites in args.sizes:
        simulate(sites, args.runs, args.workers)


if __name__ == '__main__':
    main()
