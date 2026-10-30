"""
Benchmark the corrected (spectral) noise generator to see whether a full
campaign with a theta-continuous generator is affordable.

This benchmarks the specified finite-record power-law filter and growth kernel.
DC removal, periodic boundaries and the observation duration remain part of the
protocol. Timing results do not identify a universality class or establish
agreement with an infinite-record fractional Gaussian noise covariance.

Cost model: one batched FFT of shape (n_chunk, T) per chunk of sites, so the
generation cost is O(L T log T) with numpy's vectorised pocketfft.

Usage: python code/bench_spectral.py
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bd_core import fft_noise_field, _simulate_from_field

print("=" * 96)
print("Spectral-generator benchmark  (bracketed: generation | dynamics)")
print("=" * 96)
print(f"{'L':>6} {'T':>9} {'runs':>6} {'gen s':>9} {'dyn s':>9} "
      f"{'gen s/run':>10} {'mem MB':>8}")
print("-" * 96)

for L, T, n_runs in [(64, 7680, 64), (128, 21722, 64), (256, 61440, 64),
                     (512, 173780, 32), (1024, 491520, 8)]:
    mem = L * T * 8 / 1e6
    t0 = time.time()
    try:
        eta = fft_noise_field(0.15, L, T, seed=7)
        tg = time.time() - t0
    except MemoryError:
        print(f"{L:6d} {T:9d} {n_runs:6d}   MemoryError ({mem:.0f} MB needed)")
        continue

    stride = max(1, T // 2000)
    n_w = T // stride
    check = np.unique(np.round(np.geomspace(max(16, T // 4096), T, 8)).astype(np.int64))
    prof = np.zeros((1, len(check), L), dtype=np.int32)
    ws = np.zeros((1, n_w), dtype=np.float32)
    eta3 = eta[None, :, :]
    t0 = time.time()
    _simulate_from_field(eta3, 0, L, T, stride, check, prof, ws)
    td = time.time() - t0
    print(f"{L:6d} {T:9d} {n_runs:6d} {tg:9.2f} {td:9.2f} {tg:10.3f} {mem:8.0f}")

print("-" * 96)
print("\nCompare with the FFGN path, measured on the same machine:")
print("  L=256,  T=61440,  256 runs ->   9.4 s   (0.037 s/run)")
print("  L=1024, T=655360,  96 runs ->  88.0 s   (0.92  s/run)")
print("\nIf the spectral generator is ~10x more expensive per run, a corrected")
print("campaign is still affordable at L <= 512 with ~64 runs per cell.")
