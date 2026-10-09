"""
Effective roughness exponent alpha_eff(L) at FIXED theta, and its drift with
system size.

Motivation.  The published threshold theta_c ~ 0.25 was obtained by watching the
fitted roughness exponent rise above the KPZ value 1/2.  Ales & Lopez (PRE 99,
062139 (2019)) report alpha = 0.56 +- 0.02 at theta = 0.15 for L = 2048..8192,
which is already 3 sigma above 1/2 and therefore looks like a violation of the
renormalisation-group prediction that theta < 1/6 should give alpha = 1/2.  If
that deviation is a finite-size correction, alpha_eff(theta, L) must DECREASE
towards 1/2 as L grows; if it is asymptotic, alpha_eff must be L-independent.

We compute alpha_eff from consecutive size pairs,
    alpha_eff(L_i, L_{i+1}) = log[W_sat(L_{i+1})/W_sat(L_i)] / log[L_{i+1}/L_i],
for both generators, and compare with the published numbers.

Usage: python code/alpha_drift.py
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

TH_TH = 1.0 / 6.0

ap = argparse.ArgumentParser()
ap.add_argument("--pairs", nargs="+", default=["main:FFGN", "spec:spectral"],
                help="tag:label pairs to compare, e.g. "
                     "--pairs main:FFGN spec:powerlaw fgnspec:fGn-spectrum")
args = ap.parse_args()
PAIRS = [tuple(x.split(":", 1)) for x in args.pairs]


def cells(tag, rule="bd"):
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, f"{rule}_*.npz"))):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if m:
            tok = m.group("tok")
            sign = -1.0 if tok[0] == "m" else 1.0
            out[(round(sign * int(tok[3:]) / 10000.0, 4), int(m.group("L")))] = f
    return out


def wsat(path, idx=None):
    try:
        z = np.load(path, allow_pickle=False)
        w2 = z["w_series"].astype(np.float64)
        if idx is not None:
            w2 = w2[idx]
        return float(np.sqrt(w2.mean(axis=0))[-1])
    except Exception:
        return np.nan


def alpha_eff_table(tag):
    cs = cells(tag)
    Ls = sorted({L for (_, L) in cs})
    out = {}
    for th in sorted({t for (t, L) in cs}):
        W = {}
        for L in Ls:
            if (th, L) in cs:
                w = wsat(cs[(th, L)])
                if np.isfinite(w):
                    W[L] = w
        pairs = []
        ks = sorted(W)
        for a, b in zip(ks[:-1], ks[1:]):
            pairs.append((a, b, float(np.log(W[b] / W[a]) / np.log(b / a))))
        out[th] = (W, pairs)
    return out, Ls


def alpha_bootstrap(tag, theta, n_boot=400, seed=0):
    """alpha from W_sat ~ L^alpha with a bootstrap error over independent runs.

    The naive regression error from the log-log residuals is far too optimistic
    because the run-to-run scatter of W_sat is not captured by 4-6 points.  We
    resample runs with replacement within each cell and refit.
    """
    cs = cells(tag)
    Ls = sorted({L for (t, L) in cs if t == theta})
    if len(Ls) < 3:
        return np.nan, np.nan, []
    rng = np.random.default_rng(seed)
    W0 = []
    cached = {}
    for L in Ls:
        w2 = np.load(cs[(theta, L)], allow_pickle=False)["w_series"].astype(np.float64)
        cached[L] = w2
        W0.append(float(np.sqrt(w2.mean(axis=0))[-1]))
    alpha0 = float(np.polyfit(np.log(Ls), np.log(W0), 1)[0])
    Wb = []
    for _ in range(n_boot):
        vals = []
        for L in Ls:
            w2 = cached[L]
            nr = w2.shape[0]
            idx = rng.integers(0, nr, size=nr)
            vals.append(float(np.sqrt(w2[idx].mean(axis=0))[-1]))
        Wb.append(np.polyfit(np.log(Ls), np.log(vals), 1)[0])
    return alpha0, float(np.std(Wb, ddof=1)), Ls


def main():
    print("=" * 100)
    print("alpha_eff(L) at fixed theta, from consecutive size pairs")
    print("=" * 100)

    for tag, lab in PAIRS:
        tab, Ls = alpha_eff_table(tag)
        print(f"\n--- {lab}  [tag={tag}] ---  sizes available: {Ls}")
        first = next(iter(tab.values()))[1]
        hdr = "  " + f"{'theta':>8} | " + " ".join(f"{f'{a}->{b}':>13}"
                                                   for a, b, _ in first)
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for th in sorted(tab):
            _, pairs = tab[th]
            print("  " + f"{th:8.4f} | " + " ".join(f"{v:13.4f}" for _, _, v in pairs))

    print("\n" + "=" * 100)
    print("Bootstrap error bars on the full-window alpha (W_sat ~ L^alpha)")
    print("=" * 100)
    print(f"{'generator':>10} {'theta':>8} {'alpha':>9} {'err(boot)':>10} {'L range':>14}")
    for tag, lab in PAIRS:
        cs = cells(tag)
        for th in sorted({t for (t, _) in cs}):
            if th not in (-0.2, -0.1, 0.0, 0.10, 0.1667, 0.25, 0.40):
                continue
            a, e, Ls = alpha_bootstrap(tag, th, 300)
            if np.isfinite(a):
                print(f"{lab:>10} {th:8.4f} {a:9.4f} {e:10.4f} "
                      f"{f'{min(Ls)}-{max(Ls)}':>14}")

    print("\n" + "=" * 100)
    print("Comparison with Ales & Lopez 2019 (FFGN, alpha = 0.56 +- 0.02 at theta = 0.15)")
    print("=" * 100)
    tab_ff, _ = alpha_eff_table("main")
    tab_sp, _ = alpha_eff_table("spec")
    for th_target in (0.10, 0.125, 0.15, 0.1667, 0.20, 0.25):
        best = min(tab_ff, key=lambda t: abs(t - th_target))
        if abs(best - th_target) > 0.01:
            continue
        _, pf = tab_ff[best]
        sp = tab_sp.get(best)
        f_str = ", ".join(f"{v:.4f}" for _, _, v in pf)
        s_str = ", ".join(f"{v:.4f}" for _, _, v in sp[1]) if sp else "--"
        trend = ""
        if len(pf) >= 2:
            trend = "DECREASING" if pf[-1][2] < pf[0][2] else "flat/increasing"
        print(f"\n  theta = {best:.4f}")
        print(f"    FFGN     alpha_eff by size pair: {f_str}   [{trend}]")
        print(f"    spectral alpha_eff by size pair: {s_str}")

    print("\n" + "=" * 100)
    print("Interpretation guide")
    print("=" * 100)
    print("  * theta < 1/6 and alpha_eff decreasing towards 0.5 as L grows  ->")
    print("    the published deviation is a finite-size correction and the")
    print("    renormalisation-group prediction alpha = 1/2 survives.")
    print("  * theta < 1/6 and alpha_eff L-independent and > 0.5           ->")
    print("    genuine violation of the prediction.")
    print("  * theta > 1/6 and alpha_eff decreasing                        ->")
    print("    the crossover to the long-range branch lies beyond the")
    print("    accessible sizes.")
    print("\n  CAUTION: with only 3-5 size pairs the scatter of the pair slopes is")
    print("  comparable to the deviations being tested; see the bootstrap errors")
    print("  above before drawing conclusions from any single pair.")


if __name__ == "__main__":
    main()
