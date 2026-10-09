r"""
Repair the two regressions introduced by the previous layout fix.

What went wrong
  1. Table I was given long per-column headers, e.g. "$\alpha_s$ (spectral)".
     Six such headers plus the theta column need roughly 495 pt, against 508 pt
     in a full-width table*, so the table sat right at the margin and could
     overflow.  Fix: restore a grouped header using \multicolumn + \cline in a
     PLAIN tabular -- that combination is standard LaTeX and is only dangerous
     inside revtex's \ruledtabular, which is not used here.  The table then needs
     about 300 pt and has a wide safety margin.
  2. The Appendix C protocol table was turned into a full-width table* FLOAT to
     stop it overflowing.  That fixed the width but made it drift to a page top
     and away from Appendix C, which is what "the appendix is in the wrong
     place" means.  Fix: make it a NON-floating inline centre with p{} columns
     sized in units of \columnwidth.  p{} columns wrap, so it cannot overflow,
     and being inline it stays exactly where Appendix C is.

Both repairs are checked afterwards by the tabular width estimator in
code/check_tex.py.
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


# ---------------------------------------------------------------- 1. Table I
rep(r"""\begin{tabular}{lcccccc}
\hline
$\theta$
  & $\alpha_s$ (FFGN) & $\alpha_s$ (spectral)
  & $p_{\rm flat}$ (FFGN) & $p_{\rm flat}$ (spectral)
  & $W_{\rm sat}$ (FFGN) & $W_{\rm sat}$ (spectral)\\
\hline""",
    r"""\begin{tabular}{lcccccc}
\hline
$\theta$ & \multicolumn{2}{c}{$\alpha_s$} & \multicolumn{2}{c}{$p_{\rm flat}$}
         & \multicolumn{2}{c}{$W_{\rm sat}$}\\
\cline{2-3}\cline{4-5}\cline{6-7}
 & FFGN & exact fGn & FFGN & exact fGn & FFGN & exact fGn\\
\hline""",
    "1: Table I header -> grouped, short (fits with wide margin)")

rep(r"""\theta$. $\alpha_s$ is the spectral roughness exponent, $p_{\rm flat}$ the
flat-pair fraction and $W_{\rm sat}$ the saturated width. Each quantity is given
for both generators in adjacent columns.}""",
    r"""\theta$. $\alpha_s$ is the spectral roughness exponent, $p_{\rm flat}$ the
flat-pair fraction and $W_{\rm sat}$ the saturated width. For each quantity the
left column is the FFGN (published) generator and the right column the
$\theta$-continuous exact-fGn filter.}""",
    "1: Table I caption clarified")

# ------------------------------------------------- 2. Appendix C protocol table
rep(r"""\begin{table*}[t]
\caption{\label{tab:protocol} Numerical protocol. Every quantity in the main
text is produced by these settings.}
\begin{tabular}{@{}p{0.24\textwidth}p{0.70\textwidth}@{}}""",
    r"""The settings used for every simulation reported above are listed here.
The table is set in place rather than floated, so that it stays with this
appendix.

\begin{center}
\begin{tabular}{@{}p{0.30\columnwidth}p{0.60\columnwidth}@{}}""",
    "2: protocol table -> non-floating inline centre")

rep(r"""sample matching & identical run count ($96$) at every $L$ and $\theta$, for both
  diagnostics.\\
\hline
\end{tabular}
\end{table*}""",
    r"""sample matching & identical run count ($96$) at every $L$ and $\theta$, for both
  diagnostics.\\
\hline
\end{tabular}
\end{center}""",
    "2: protocol table footer -> inline")

open(p, "w", encoding="utf-8", newline="\n").write(src)
print(f"\nmanuscript/main.tex: {n0} -> {len(src)} chars")
print(f"table* count: {src.count(chr(92)+'begin{table*}')}, "
      f"table count: {src.count(chr(92)+'begin{table}')}")
