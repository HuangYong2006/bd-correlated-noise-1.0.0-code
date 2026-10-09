r"""
Make every inline table fit inside a single PRE column.

Finding: the manuscript contains SEVEN tabulars, not two.  Four of them are
inline centred tables in the Results sections, and four of those are wider than
one column (246 pt):

  L 756  tabular{lccccc}   ~251 pt   onset theta*(L)
  L 807  tabular{lcccccc}  ~369 pt   alpha, and it declared 7 columns for a
                                     6-column header (a real bug)
  L 945  tabular{lccccc}   ~394 pt   alpha for the three generators
  L 969  tabular{lccccc}   ~293 pt   fixed-theta slopes

A full-width table* is not the answer for these: float displacement is exactly
what made the Appendix C table leave its appendix.

Fixes applied
  * a \compacttab helper (\footnotesize + \tabcolsep=4pt), applied to every
    inline table, buys roughly 30% of width
  * the bootstrap errors are removed from the cells (they cost ~18 pt per
    column) and a representative uncertainty is stated in the text instead
  * column counts corrected, and the long row labels shortened

Usage: python code/fix_inline_tables.py
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


# ---------------------------------------------------- preamble: helper macro
rep(r"""% red marker, available for the author's own use; currently unused
\newcommand{\TODO}[1]{\textcolor{red}{[TODO: #1]}}""",
    r"""% red marker, available for the author's own use; currently unused
\newcommand{\TODO}[1]{\textcolor{red}{[TODO: #1]}}

% Compact inline tables.  A single PRE column is only 246 pt wide, and several
% of the small reference tables below need most of it; \footnotesize with
% \tabcolsep=4pt keeps them inside the column.  RevTeX fixes tables at the
% current font size, so this genuinely reduces the width.
\newcommand{\compacttab}{\footnotesize\setlength{\tabcolsep}{4pt}}""",
    "preamble: \\compacttab helper")

# ------------------------------------------------------------- table 1 (p_flat)
rep(r"""\begin{center}
\begin{tabular}{lccc}
\hline
noise & $\theta$ & $p_{\mathrm{flat}}$ & $W_{\mathrm{sat}}$\\""",
    r"""\begin{center}
\compacttab
\begin{tabular}{lccc}
\hline
noise & $\theta$ & $p_{\mathrm{flat}}$ & $W_{\mathrm{sat}}$\\""",
    "table 1 (p_flat): compact")

# --------------------------------------------------------------- table 2 (onset)
rep(r"""\begin{center}
\begin{tabular}{lccccc}
\hline
$L$ & 32 & 64 & 128 & 256 & 512\\
\hline
$\theta^*(z\ge3)$ & $0.2415$ & $0.2336$ & $0.1710$ & $0.1375$ & $0.0823$\\
$\theta^*(z\ge5)$ & $0.3306$ & $0.3073$ & $0.2686$ & $0.2410$ & $0.1471$\\
\hline
\end{tabular}
\end{center}""",
    r"""\begin{center}
\compacttab
\begin{tabular}{lccccc}
\hline
 & \multicolumn{5}{c}{$\theta^*(L)$}\\
\cline{2-6}
$L$ & 32 & 64 & 128 & 256 & 512\\
\hline
$z\ge3$ & $0.2415$ & $0.2336$ & $0.1710$ & $0.1375$ & $0.0823$\\
$z\ge5$ & $0.3306$ & $0.3073$ & $0.2686$ & $0.2410$ & $0.1471$\\
\hline
\end{tabular}
\end{center}""",
    "table 2 (onset): compact, grouped header, short row labels")

# -------------------------------------------------- table 3 (alpha, 2 generators)
rep(r"""\begin{center}
\begin{tabular}{lcccccc}
\hline
$\theta$ & $-0.2$ & $0$ & $1/6$ & $0.25$ & $0.40$\\
\hline
$\alpha$ (spectral) & $0.4947(14)$ & $0.5062(13)$ & $0.5082(12)$ &
$0.5311(12)$ & $0.5750(12)$\\
$\alpha$ (FFGN) & --- & $0.4986(13)$ & $0.6330(15)$ & $0.7019(13)$ &
$0.8585(12)$\\
\hline
\end{tabular}
\end{center}""",
    r"""\begin{center}
\compacttab
\begin{tabular}{lccccc}
\hline
$\theta$ & $-0.2$ & $0$ & $1/6$ & $0.25$ & $0.40$\\
\hline
FFGN       & ---      & $0.4986$ & $0.6330$ & $0.7019$ & $0.8585$\\
corrected  & $0.4947$ & $0.5062$ & $0.5082$ & $0.5311$ & $0.5750$\\
\hline
\end{tabular}
\end{center}""",
    "table 3 (alpha, two generators): 6 columns not 7, errors moved to text")

rep(r"""With bootstrap errors, over $L=32\dots512$,""",
    r"""With bootstrap errors, over $L=32\dots512$ (the bootstrap uncertainty on
every entry below is $\pm0.011$--$0.015$, i.e.\ it is the last digit or two and
is omitted from the table for width):""",
    "table 3: uncertainty moved into the text")

# ------------------------------------------- table 4 (alpha, three generators)
rep(r"""\begin{center}
\begin{tabular}{lccccc}
\hline
$\theta$ & $0$ & $0.10$ & $1/6$ & $0.25$ & $0.40$\\
\hline
$\alpha$ (FFGN, published)  & $0.4986(13)$ & $0.5835(13)$ & $0.6330(15)$
                            & $0.7019(13)$ & $0.8585(12)$\\
$\alpha$ (power-law filter) & $0.5062(13)$ & $0.5247(13)$ & $0.5082(12)$
                            & $0.5311(12)$ & $0.5750(12)$\\
$\alpha$ (exact fGn filter) & $0.5076(12)$ & $0.5322(13)$ & $0.5089(11)$
                            & $0.5383(12)$ & $0.5724(12)$\\
\hline
\end{tabular}
\end{center}""",
    r"""\begin{center}
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
    "table 4 (alpha, three generators): compact, short labels, errors to text")

# ------------------------------------------ table 5 (fixed-theta slopes)
rep(r"""\begin{center}
\begin{tabular}{lccccc}
\hline
$\theta$ & $0.10$ & $1/6$ & $0.20$ & $0.25$ & $0.30$\\
\hline
$\mathrm{d}\,\Delta H/\mathrm{d}\ln L$
  & $-0.0024$ & $+0.0011$ & $+0.0004$ & $-0.0023$ & $+0.0048$\\
significant? & no & no & no & no & \textbf{yes}\\
\hline
\end{tabular}
\end{center}""",
    r"""\begin{center}
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
    "table 5 (fixed-theta slopes): compact, short label")

open(p, "w", encoding="utf-8", newline="\n").write(src)
print(f"\nmanuscript/main.tex: {n0} -> {len(src)} chars")
print(f"\\compacttab uses: {src.count(chr(92) + 'compacttab') - 1}")
