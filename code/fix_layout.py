r"""
Fix the layout defects found in the compiled PDF (manuscript/overleaf_bd_paper.pdf).

Diagnosis from code/overfull.py -- 12 overfull lines, two causes:

  A. Appendix C's protocol table used tabular{ll}: the second column is long
     prose in an unbreakable 'l' column, so it ran 30-54 pt past the right
     margin and off the page edge (x up to 616 pt against a 612 pt page).
  B. Equation (eq:ffgncorr) carried "H = theta + 1/2" on the same line and
     overflowed its column by 29.5 pt.

Also fixed here:
  C. Table I is promoted to a full-width table* with explicit column labels and
     a plain tabular (no ruledtabular, no multicolumn/cline), which is both
     journal-appropriate and free of the "Misplaced \noalign" risk.
  D/E. The two remaining red \TODO markers are replaced by ordinary prose and a
     plain placeholder, so the compiled PDF contains no red placeholder text.

Only manuscript/main.tex is edited; the Overleaf package is rebuilt afterwards.

Usage: python code/fix_layout.py
"""

from __future__ import annotations

import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
p = os.path.join(ROOT, "manuscript", "main.tex")
src = open(p, encoding="utf-8").read()
n0 = len(src)
report = []


def rep(old, new, label):
    global src
    c = src.count(old)
    if c != 1:
        raise SystemExit(f"anchor count {c} (need 1) for {label}")
    src = src.replace(old, new)
    report.append(f"  {label}")


# ---------------------------------------------------------------- B: equation
rep(r"""\begin{equation}
  C_{\mathrm{FFGN}}(\tau) \;\sim\; A_{\mathrm{FFGN}}\,\tau^{2H-2},
  \qquad H = \theta+\tfrac12 ,
  \label{eq:ffgncorr}
\end{equation}
with the coefficient obtained from the geometric-grid integral""",
    r"""\begin{equation}
  C_{\mathrm{FFGN}}(\tau) \;\sim\; A_{\mathrm{FFGN}}\,\tau^{-2\kappa},
  \qquad \kappa = \tfrac12-\theta ,
  \label{eq:ffgncorr}
\end{equation}
with the coefficient obtained from the geometric-grid integral""",
    "B: eq:ffgncorr shortened (kappa notation, fits the column)")

# the exponent 2H-2 is used elsewhere; keep consistency by stating H there
rep(r"""with a \emph{single} white sequence $\zeta_t$ driving all components. Because
Eq.~\eqref{eq:ffgn} is a linear filter of one white sequence, its correlation
function is""",
    r"""with a \emph{single} white sequence $\zeta_t$ driving all components. Because
Eq.~\eqref{eq:ffgn} is a linear filter of one white sequence, its correlation
function decays as a power law with a $\theta$-dependent amplitude,""",
    "B: lead-in sentence adjusted")

# ------------------------------------------------------- C: Table I -> table*
rep(r"""\begin{table}[h]
\caption{\label{tab:morph} Exponents and morphology at $L=256$ for the FFGN
(published) and spectral ($\theta$-continuous) generators at matched nominal
$\theta$. $\alpha_s$ is the spectral roughness exponent, $p_{\rm flat}$ the
flat-pair fraction and $W_{\rm sat}$ the saturated width.}
\begin{ruledtabular}
\begin{tabular}{lcccccc}
$\theta$ & \multicolumn{2}{c}{$\alpha_s$} & \multicolumn{2}{c}{$p_{\rm flat}$}
         & \multicolumn{2}{c}{$W_{\rm sat}$}\\
\cline{2-3}\cline{4-5}\cline{6-7}
 & FFGN & spec & FFGN & spec & FFGN & spec\\
\hline""",
    r"""\begin{table*}[t]
\caption{\label{tab:morph} Exponents and morphology at $L=256$ for the FFGN
(published) and spectral ($\theta$-continuous) generators at matched nominal
$\theta$. $\alpha_s$ is the spectral roughness exponent, $p_{\rm flat}$ the
flat-pair fraction and $W_{\rm sat}$ the saturated width. Each quantity is given
for both generators in adjacent columns.}
\begin{tabular}{lcccccc}
\hline
$\theta$
  & $\alpha_s$ (FFGN) & $\alpha_s$ (spectral)
  & $p_{\rm flat}$ (FFGN) & $p_{\rm flat}$ (spectral)
  & $W_{\rm sat}$ (FFGN) & $W_{\rm sat}$ (spectral)\\
\hline""",
    "C: Table I header -> full-width table*, explicit labels")

rep(r"""\end{tabular}
\end{ruledtabular}
\end{table}""",
    r"""\hline
\end{tabular}
\end{table*}""",
    "C: Table I footer")

# ------------------------------------------------ A: appendix protocol table
rep(r"""\begin{center}
\begin{tabular}{ll}
\hline
quantity & value\\
\hline
FFGN track (published generator) & $L=64\dots512$, $T=15L^{3/2}$, up to $256$ runs\\
spectral track (corrected) & $L=32\dots256$ ($512$ for selected $\theta$), $T=15L^{3/2}$,
  $192$ runs ($96$ at $L=512$)\\
$\theta$ grid (FFGN) & $0,0.0025,\dots,0.4$ (22 values)\\
$\theta$ grid (spectral) & $0$, $\pm$ negative controls, and $0.01$--$0.02$ steps
  across $[0.10,0.225]$ (26 values)\\
snapshots per run & $32$ geometrically spaced times\\
width series & $2\times10^4$ equidistant samples\\
saturation check & $W$ flat over the last decade of $t$\\
Spectral FFT length & 5-smooth (rounded up from $15L^{3/2}$), threaded \texttt{scipy.fft}\\
Hellinger null floor & split-half resampling, $60$--$80$ splits\\
Isolation Forest & $200$--$400$ trees, FPR $\alpha=0.02$ calibrated on held-out
  $\theta=0$ data, $3$--$16$ train/calibration resplits\\
iForest sample matching & identical run count ($96$) at every $L$ and $\theta$\\
\hline
\end{tabular}
\end{center}""",
    r"""\begin{table*}[t]
\caption{\label{tab:protocol} Numerical protocol. Every quantity in the main
text is produced by these settings.}
\begin{tabular}{@{}p{0.24\textwidth}p{0.70\textwidth}@{}}
\hline
quantity & value\\
\hline
FFGN track (published generator) & $L=64\dots512$, $T=15L^{3/2}$ sweeps, up to
  $256$ independent realisations per $(\theta,L)$ cell.\\
spectral track (corrected generator) & $L=32\dots256$ over the full $\theta$
  grid ($512$ for a subset of $\theta$), $T=15L^{3/2}$ sweeps, $192$
  realisations per cell ($96$ at $L=512$).\\
$\theta$ grid (FFGN) & $0,\,0.0025,\,0.005,\dots,0.40$, $22$ values.\\
$\theta$ grid (spectral) & $0$, negative controls at $-0.05,-0.1,-0.2$, and
  $0.010$--$0.020$ steps across $[0.10,0.225]$; $26$ values in total.\\
snapshots per run & $32$ geometrically spaced times.\\
width series & $2\times10^{4}$ equidistant samples.\\
saturation check & $W$ flat over the last decade of $t$.\\
spectral FFT length & 5-smooth (rounded up from $15L^{3/2}$), threaded
  \texttt{scipy.fft}.\\
Hellinger null floor & split-half resampling over runs, $60$--$80$ splits.\\
Isolation Forest & $200$--$400$ trees, false-positive rate $\alpha=0.02$
  calibrated on held-out $\theta=0$ data, $3$--$16$ train/calibration
  resplits.\\
sample matching & identical run count ($96$) at every $L$ and $\theta$, for both
  diagnostics.\\
\hline
\end{tabular}
\end{table*}""",
    "A: protocol table -> full-width table* with wrapping p{} columns")

# ------------------------------------------------------------------ D: TODO 1
rep(r"""\TODO{verify the novelty-detection algorithm used in Ref.~\cite{WangXia2025}
once the full text is accessible, and adjust this paragraph accordingly.}""",
    r"""One limitation should be stated at the outset. Reference~\cite{WangXia2025}
is paywalled and we were unable to establish which novelty-detection algorithm
it employs, so the methodological critique developed below applies to the onset
construction itself rather than to that particular choice of detector.""",
    "D: novelty-detection TODO -> prose")

# ------------------------------------------------------------------ E: TODO 2
rep(r"""\begin{acknowledgments}
\TODO{acknowledgments}
\end{acknowledgments}""",
    r"""\begin{acknowledgments}
% FILL IN: funding, computing resources, and colleagues to thank.
The authors thank \dots\ for discussions. Computations were performed on
\dots.
\end{acknowledgments}""",
    "E: acknowledgments TODO -> plain placeholder")

open(p, "w", encoding="utf-8", newline="\n").write(src)
print("manuscript/main.tex patched:")
print("\n".join(report))
print(f"\n{n0} -> {len(src)} chars")
print(f"remaining \\TODO uses: {src.count(chr(92) + 'TODO{')}")
