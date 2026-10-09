r"""
Correction to fix_layout.py.

The first pass replaced (2H-2, H=theta+1/2) by (-2*kappa, kappa=1/2-theta) in
eq:ffgncorr.  That removed the only place where H is defined, while H is still
used by eq:affgn (B^{1-H}, Gamma(2-2H), Gamma(3-2H)) and elsewhere.  Restore H
but keep the equation short enough for the column by moving the definition of H
into the preceding sentence.

Usage: python code/fix_layout2.py
"""

from __future__ import annotations

import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
p = os.path.join(ROOT, "manuscript", "main.tex")
src = open(p, encoding="utf-8").read()
n0 = len(src)


def rep(old, new, label):
    global src
    c = src.count(old)
    if c != 1:
        raise SystemExit(f"anchor count {c} (need 1) for {label}")
    src = src.replace(old, new)
    print(f"  {label}")


rep(r"""function decays as a power law with a $\theta$-dependent amplitude,""",
    r"""function decays as a power law, with $H=\theta+\tfrac12$ and a
$\theta$-dependent amplitude,""",
    "H restored in the lead-in sentence")

rep(r"""  C_{\mathrm{FFGN}}(\tau) \;\sim\; A_{\mathrm{FFGN}}\,\tau^{-2\kappa},
  \qquad \kappa = \tfrac12-\theta ,
  \label{eq:ffgncorr}""",
    r"""  C_{\mathrm{FFGN}}(\tau) \;\sim\; A_{\mathrm{FFGN}}\,\tau^{2H-2},
  \label{eq:ffgncorr}""",
    "eq:ffgncorr back to the 2H-2 form, short enough for the column")

open(p, "w", encoding="utf-8", newline="\n").write(src)
print(f"\n{n0} -> {len(src)} chars")

# report every remaining \TODO occurrence with its line, so nothing is missed
print("\nremaining \\TODO occurrences (should all be inside % comments):")
for i, ln in enumerate(src.split("\n"), 1):
    if "\\TODO" in ln:
        stripped = ln.strip()
        where = "COMMENT" if stripped.startswith("%") else "!!! IN BODY !!!"
        print(f"  L{i:5d} [{where}] {stripped[:88]}")
