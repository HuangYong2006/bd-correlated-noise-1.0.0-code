r"""
Transpose the two remaining tight tables.

Both are 6-column tables at 88-90% of a single PRE column, which is too close
to the edge given that the width estimator carries a +/-20% error.  Transposing
them puts the generators (or the quantities) in columns and the parameter in
rows, which both fits comfortably and reads better, since the comparison of
interest -- between generators -- becomes a row-wise comparison.

Usage: python code/transpose_tables.py
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


# ------------------------------- alpha table: theta becomes the row index -----
rep(r"""\begin{center}
\compacttab
\begin{tabular}{lccccc}
\hline
$\theta$ & $0$ & $1/6$ & $0.20$ & $0.25$ & $0.40$\\
\hline
FFGN (published) & $0.4986$ & $0.6330$ & $0.6855$ & $0.7019$ & $0.8585$\\
power-law filter & $0.5062$ & $0.5082$ & $0.5269$ & $0.5311$ & $0.5750$\\
exact fGn filter & $0.5076$ & $0.5089$ & $0.5368$ & $0.5383$ & $0.5724$\\
\hline
\end{tabular}
\end{center}""",
                                    
    "alpha table (3 generators): transposed to 4 columns")

# --------------------------- fixed-theta slope table: one row per theta -------
rep(r"""\begin{center}
\compacttab
\begin{tabular}{lccccc}
\hline
$\theta$ & $0.10$ & $1/6$ & $0.20$ & $0.25$ & $0.30$\\
\hline
slope & $-0.0024$ & $+0.0011$ & $+0.0004$ & $-0.0023$ & $+0.0048$\\
significant? & no & no & no & no & \textbf{yes}\\
\hline
\end{tabular}
\end{center}""",
    r"""\begin{center}
\compacttab
\begin{tabular}{lcc}
\hline
$\theta$ & slope & sig.?\\
\hline
$0.10$ & $-0.0024$ & no\\
$1/6$  & $+0.0011$ & no\\
$0.20$ & $+0.0004$ & no\\
$0.25$ & $-0.0023$ & no\\
$0.30$ & $+0.0048$ & \textbf{yes}\\
\hline
\end{tabular}
\end{center}""",
    "fixed-theta slope table: transposed to 3 columns")

open(p, "w", encoding="utf-8", newline="\n").write(src)
print(f"\nmanuscript/main.tex: {n0} -> {len(src)} chars")
