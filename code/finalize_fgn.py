r"""
Update the fig:drift caption for the three-generator version, verify the
manuscript structurally, and rebuild the Overleaf package.
"""

from __future__ import annotations

import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.join(ROOT, "code")
p = os.path.join(ROOT, "manuscript", "main.tex")
src = open(p, encoding="utf-8").read()
n0 = len(src)

old_cap = r"""  \caption{\label{fig:drift} The decisive tests. (a) $\alpha_{\rm eff}$ from
  consecutive size pairs versus $L$, for $\theta\in\{0.10,1/6,0.20,0.25\}$
  (dark to light), solid for the $\theta$-continuous generator and dashed for
  FFGN. The corrected $\alpha_{\rm eff}$ falls towards $1/2$ as $L$ grows for
  $\theta\ge1/6$; the FFGN one does not. (b) The fixed-$\theta$ test:
  $\Delta H(\theta,L)$ versus $L$ for the same values of $\theta$. A genuine
  crossover requires the deviation to \emph{grow} with $L$; it is instead
  $L$-independent within errors (symbols as in (a)). (c) $\alpha(\theta)$ with
  bootstrap errors obtained by resampling runs, against the Medina one-loop and
  Squizzato--Canet long-range branches. (d) The split-half null floor of
  $\Delta H$ versus $L$: it falls by a factor of $4.5$ between $L=32$ and $512$,
  which is what makes a fixed-significance onset drift even when the observable
  is $L$-independent.}"""

new_cap = r"""  \caption{\label{fig:drift} The decisive tests, for the two contested values
  $\theta=0.20$ (circles) and $\theta=0.25$ (squares), comparing three noise
  generators: FFGN (dashed), the power-law spectral filter (solid) and the
  exact-fractional-Gaussian-noise spectral filter (dash-dotted).
  (a) $\alpha_{\rm eff}$ from consecutive size pairs versus $L$. Both corrected
  filters fall towards $1/2$ as $L$ grows; the FFGN value does not.
  (b) The fixed-$\theta$ test: $\Delta H(\theta,L)$ versus $L$. A genuine
  crossover requires the deviation to \emph{grow} with $L$; at these two values
  of $\theta$ it is $L$-independent within errors for both corrected filters.
  (c) $\alpha(\theta)$ with bootstrap errors obtained by resampling runs,
  against the Medina one-loop and Squizzato--Canet long-range branches. The two
  corrected filters agree to within $0.008$ everywhere, while FFGN lies
  $0.16$--$0.29$ higher. (d) The split-half null floor of $\Delta H$ versus $L$:
  it falls by a factor of $4.5$ between $L=32$ and $512$, which is what makes a
  fixed-significance onset drift even when the observable is $L$-independent.}"""

if src.count(old_cap) != 1:
    raise SystemExit(f"caption anchor count {src.count(old_cap)} (need 1)")
src = src.replace(old_cap, new_cap)
open(p, "w", encoding="utf-8", newline="\n").write(src)
print(f"caption updated; main.tex {n0} -> {len(src)} chars")

print("\n--- structural check: source ---")
r = subprocess.run([sys.executable, os.path.join(HERE, "check_tex.py"),
                    os.path.join(ROOT, "manuscript", "main.tex")],
                   capture_output=True, text=True, encoding="utf-8")
print("\n".join(r.stdout.splitlines()[:11]))

print("\n--- rebuild overleaf package ---")
r = subprocess.run([sys.executable, os.path.join(HERE, "build_overleaf.py")],
                   capture_output=True, text=True, encoding="utf-8")
print(r.stdout.strip() or r.stderr.strip())

print("\n--- structural check: overleaf ---")
r = subprocess.run([sys.executable, os.path.join(HERE, "check_tex.py"),
                    os.path.join(ROOT, "manuscript", "overleaf", "main.tex")],
                   capture_output=True, text=True, encoding="utf-8")
print("\n".join(r.stdout.splitlines()[:11]))
