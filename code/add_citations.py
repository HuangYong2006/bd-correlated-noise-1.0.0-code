r"""
Add the seven bibliography entries that were present in refs.bib but never
cited, at the points where they belong.  Applied to manuscript/main.tex, after
which the Overleaf package is rebuilt.

  FamilyVicsek1985 + BarabasiStanley1995 -> the scaling Ansatz / exponent
      definitions in the simulation-protocol subsection
  Hellinger1909 + Bhattacharyya1943      -> the Hellinger-distance definition
  HuHaoXia2023                           -> the numerical literature review
      (they report alpha_loc ~ alpha < alpha_s for large theta and state that
      continuum and discrete growth systems need not share a theta)
  ChangXia2025                           -> spatiotemporal correlations
  Szendro2007                            -> the columnar-disorder theta = 1/2 limit
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
        raise SystemExit(f"anchor not unique ({c}) for {label}: {old[:70]!r}")
    src = src.replace(old, new)
    print(f"  {label:44s} ok")


# 1. scaling Ansatz / exponents -------------------------------------------
rep(r"""The KPZ saturation time of model \eqref{eq:bd} was measured to be""",
    r"""The interface width follows the standard Family--Vicsek scaling Ansatz
\cite{FamilyVicsek1985,BarabasiStanley1995},
$W(L,t)=L^{\alpha}f(t/L^{z})$ with $z=\alpha/\beta$, so that
$W\sim t^{\beta}$ in the growth regime and $W_{\mathrm{sat}}\sim L^{\alpha}$ at
saturation; $\alpha$, $\beta$ and $z$ are the quantities whose $\theta$-dependence
is at issue here.

The KPZ saturation time of model \eqref{eq:bd} was measured to be""",
    "Family-Vicsek scaling Ansatz")

# 2. Hellinger / Bhattacharyya --------------------------------------------
rep(r"""the Bhattacharyya coefficient $\sum_k\sqrt{p_kq_k}$ being the fidelity between
$P$ and $Q$.""",
    r"""the Bhattacharyya coefficient $\sum_k\sqrt{p_kq_k}$ being the fidelity between
$P$ and $Q$ \cite{Hellinger1909,Bhattacharyya1943}.""",
    "Hellinger / Bhattacharyya citations")

# 3. numerical literature review ------------------------------------------
rep(r"""Hellinger distance and a novelty detector. The non-perturbative renormalisation
group authors describe the situation frankly:""",
    r"""Hellinger distance and a novelty detector. Hu, Hao and Xia
\cite{HuHaoXia2023} obtain $\alpha_{\mathrm{loc}}\simeq\alpha<\alpha_s$ for
large $\theta$ and conclude that continuum and discrete growth systems ``do not
belong to the same universality class with the same temporal correlation
exponent''; Chang and Xia \cite{ChangXia2025} extend the distribution analysis
to jointly spatial and temporal correlations; and the strong-correlation limit
$\theta\to1/2$ has been treated as columnar disorder \cite{Szendro2007}. The
non-perturbative renormalisation group authors describe the situation frankly:""",
    "HuHaoXia2023 / ChangXia2025 / Szendro2007")

open(p, "w", encoding="utf-8", newline="\n").write(src)
print(f"\nmain.tex: {n0} -> {len(src)} chars")
