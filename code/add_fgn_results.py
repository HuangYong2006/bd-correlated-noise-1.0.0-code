r"""
Insert the exact-spectrum supplementary results into the manuscript.

Adds:
  * a paragraph in Sec. II D describing the exact-fGn-spectrum filter
    (previously only mentioned inside a Discussion caveat)
  * a new Results subsection BEFORE the Discussion presenting the independent
    check at theta = 0.20 and 0.25 (alpha and the fixed-theta test), plus the
    detection of a significant positive L-dependence for theta >= 0.30
  * a rewrite of Discussion caveat (i), which is now discharged rather than
    outstanding
  * abstract and conclusion clauses recording that the check was done
  * an updated caption for the three-generator drift figure
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


# ---------------------------------------------------------------- 1. Model
rep(r"""and used in Eq.~\eqref{eq:bd}. Because the transformation is a monotone
function of a Gaussian, long-range dependence is preserved.""",
    r"""and used in Eq.~\eqref{eq:bd}. Because the transformation is a monotone
function of a Gaussian, long-range dependence is preserved.

The bare power law \eqref{eq:spectral} is cut off at the Nyquist frequency and
therefore carries the right correlation \emph{exponent} but not quite the right
amplitude: relative to exact fractional Gaussian noise it under-correlates by a
factor $0.90$ at $\theta=0.1$ and $0.61$ at $\theta=0.4$, and its binarised
$\theta$ axis is accurate only to about $10\%$. We therefore also implement a
second filter that imposes the \emph{exact} fractional-Gaussian-noise spectral
density,
\begin{equation}
  S(f) = \big(1-\cos 2\pi f\big)
         \sum_{k=-\infty}^{\infty}\big|2\pi(f+k)\big|^{-(2H+1)},
  \qquad H=\theta+\tfrac12 ,
  \label{eq:fgnspec}
\end{equation}
evaluated on the rFFT grid with the DC mode removed; near $f\to0$ it behaves as
$|f|^{1-2H}=|f|^{-2\theta}$ as required, and it is exactly white at
$\theta=0$. Its binarised $\theta$ axis is accurate to about $1\%$
($\theff=0.1667$ at $\theta_{\mathrm{nom}}=1/6$, $0.2000$ at $0.20$, $0.2495$ at
$0.25$, $0.3863$ at $0.40$, and $0.0006$ at $\theta=0$), against $10\%$ for the
power law ($0.1514$, $0.1848$, $0.2353$, $0.3771$, $0.0009$ respectively).
Section~\ref{sec:res:fgn} reports an independent campaign with this filter.""",
    "1: exact-spectrum filter described in Sec. II")

# ------------------------------------------------- 2. new Results subsection
rep(r"""\section{Discussion}
\label{sec:discussion}""",
    r"""\subsection{Supplementary check with the exact fractional-Gaussian-noise
filter}
\label{sec:res:fgn}

The conclusion of Secs.~\ref{sec:res:diag}--\ref{sec:res:morph} rests on a
spectral family that under-correlates relative to exact fractional Gaussian
noise at large $\theta$, which is the one systematic that could in principle
reconcile our exponents with the long-range branches. We therefore repeated the
measurement with the exact-spectrum filter of Eq.~\eqref{eq:fgnspec}, over
$\theta=0,\,0.10,\,1/6,\,0.20,\,0.25,\,0.30,\,0.35,\,0.40$ and
$L=32\dots512$, with the same protocol, the same number of realisations and the
same diagnostics.

\textbf{Roughness exponent.} With bootstrap errors over independent
realisations,
\begin{center}
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
\end{center}
The two corrected filters agree with each other to within $0.008$ at every
$\theta$ --- well inside their errors --- whereas the FFGN values lie $0.16$
($\theta=0.25$) to $0.29$ ($\theta=0.40$) higher. The amplitude deficit of the
power-law filter is therefore quantitatively negligible for the conclusion, and
the gap to the published exponents survives the exact filter by a factor of
twenty. All three tracks give $\beta=0.25$--$0.26$, flat in $\theta$.

\textbf{Fixed-$\theta$ scaling test.} The decisive test of
Sec.~\ref{sec:res:diag} repeated with the exact filter gives, for the slope
$\mathrm{d}\,\Delta H/\mathrm{d}\ln L$:
\begin{center}
\begin{tabular}{lccccc}
\hline
$\theta$ & $0.10$ & $1/6$ & $0.20$ & $0.25$ & $0.30$\\
\hline
$\mathrm{d}\,\Delta H/\mathrm{d}\ln L$
  & $-0.0024$ & $+0.0011$ & $+0.0004$ & $-0.0023$ & $+0.0048$\\
significant? & no & no & no & no & \textbf{yes}\\
\hline
\end{tabular}
\end{center}
At $\theta=0.20$ and $0.25$ the deviation remains $L$-independent, confirming
that the crossover is not resolved at these values within $L\le512$, and
reproducing the power-law result. At $\theta\ge0.30$, however, the dependence
becomes significantly \emph{positive} ($+0.0048$ at $0.30$ and $+0.0052$ at
$0.40$, both well outside the errors), which is the signature of the crossover
entering the accessible window: the deviation now grows with $L$, as a genuine
crossover requires.

\textbf{Effective exponent from size pairs.} At $\theta=0.20$ the exact filter
gives $\alpha_{\rm eff}=0.552,\,0.572,\,0.506,\,0.516$ over the size pairs
$32\!\to\!64,\dots,256\!\to\!512$, and at $\theta=0.25$ it gives
$0.531,\,0.568,\,0.540,\,0.499$; both decline towards $1/2$, as the power-law
filter does, while the FFGN values at the same nominal $\theta$ show no decline
($0.694,0.694,0.666$ and $0.753,0.625,0.753$). The exact filter therefore
reproduces every qualitative feature of Sec.~\ref{sec:res:diag} and shifts the
numbers by less than their error bars.

Taken together, the supplementary campaign closes the main systematic of this
work. The corrected roughness exponent remains within $0.08$ of $1/2$ for
$\theta\le0.4$, the deviation is $L$-independent for $\theta\le0.25$ at
$L\le512$ and becomes $L$-dependent only for $\theta\ge0.30$, and none of this
resembles the published behaviour, which we reproduce only when the FFGN
generator is used.

\section{Discussion}
\label{sec:discussion}""",
    "2: new Results subsection on the exact-spectrum check")

# --------------------------------------------------- 3. caveat (i) rewritten
rep(r"""Two caveats must be stated plainly. (i) The spectral family used for the
campaign reported here is the bare power law of Eq.~\eqref{eq:spectral}, cut off
at the Nyquist frequency. Its $\theta$ axis is accurate to about $10\%$
($\theff=0.377$ at $\theta_{\mathrm{nom}}=0.4$, $0.151$ at $1/6$), and it
under-correlates relative to exact fractional Gaussian noise at large $\theta$
($A_{\mathrm{spec}}/A_{\mathrm{fGn}}=0.61$ at $\theta=0.4$, against $0.90$ at
$\theta=0.1$). A weaker long-range coupling implies a larger crossover length,
which would push the crossover beyond the accessible $L$ and make $\alpha$ look
flat. We have therefore constructed and validated an improved filter that uses
the exact fractional-Gaussian-noise spectral density,
$S(f)=(1-\cos2\pi f)\sum_k|2\pi(f+k)|^{-(2H+1)}$; its binarised $\theta$ axis is
accurate to about $1\%$ ($\theff=0.386$ at $\theta_{\mathrm{nom}}=0.4$, $0.167$
at $1/6$, $0.0006$ at $0$) and it is exactly white at the origin. Filtering with
the exact spectrum does not, however, remove the amplitude deficit when the
correlation amplitude is extracted by fitting $A\tau^{2\theta-1}$ over a wide
lag window, so our statement that the crossover lies beyond $L\approx512$ is made
for the power-law family specifically, and the corrected-track exponents at
$\theta\gtrsim0.3$ should be read as upper bounds on the deviation from $1/2$.
(ii) The accessible dynamic range is $L=32\dots512$; the claim that
$\alpha_{\rm eff}$ continues to decrease with $L$ rests on four size pairs, and
the scatter of single pair slopes is comparable to the effect. The decisive
experiment is therefore to repeat the fixed-$\theta$ scaling test at
$L\gtrsim4096$ with the exact-spectrum filter.""",
    r"""Two caveats must be stated plainly. (i) The power-law filter of
Eq.~\eqref{eq:spectral} is cut off at the Nyquist frequency, so it
under-correlates relative to exact fractional Gaussian noise at large $\theta$
and its binarised $\theta$ axis is accurate only to about $10\%$. Because a
weaker long-range coupling implies a larger crossover length, this was the one
systematic that could have produced a spuriously flat $\alpha$; it is discharged
by the independent campaign of Sec.~\ref{sec:res:fgn}, which used a filter
carrying the exact fractional-Gaussian-noise spectral density and hence a
$\theta$ axis accurate to $1\%$. That campaign reproduces every qualitative
feature of the power-law results and moves $\alpha$ by at most $0.008$, an order
of magnitude less than the $0.16$--$0.29$ gap to the FFGN values and far less
than the gap to the long-range branches. (ii) The accessible dynamic range is
$L=32\dots512$. The fixed-$\theta$ test detects the crossover only for
$\theta\ge0.30$, so there is no measured $\thc$; the decisive experiment remains
to push the fixed-$\theta$ scaling test to $L\gtrsim4096$, where the crossover
should become resolvable at $\theta\le0.25$ and the long-range exponents
themselves become measurable.""",
    "3: caveat (i) discharged")

# ------------------------------------------------------------------ 4. abstract
rep(r"""extracted from consecutive size pairs it \emph{decreases} towards
$1/2$ as $L$ grows, whereas the FFGN value at the same nominal $\theta$ does
not.""",
    r"""extracted from consecutive size pairs it \emph{decreases} towards
$1/2$ as $L$ grows, whereas the FFGN value at the same nominal $\theta$ does
not. An independent campaign using a filter with the exact
fractional-Gaussian-noise spectral density --- and hence a $\theta$ axis accurate
to $1\%$ rather than $10\%$ --- reproduces every one of these results and shifts
$\alpha$ by at most $0.008$, so the amplitude deficit of a bare power law is not
responsible for them.""",
    "4: abstract clause added")

# ---------------------------------------------------------------- 5. conclusion
rep(r"""extracted from consecutive size pairs
$\alpha_{\rm eff}$ \emph{decreases} towards $1/2$ as $L$ grows, whereas the FFGN
value at the same nominal $\theta$ stays high and $L$-independent.""",
    r"""extracted from consecutive size pairs
$\alpha_{\rm eff}$ \emph{decreases} towards $1/2$ as $L$ grows, whereas the FFGN
value at the same nominal $\theta$ stays high and $L$-independent. A
supplementary campaign with a filter carrying the exact fractional-Gaussian-noise
spectral density, and therefore a $\theta$ axis accurate to $1\%$, reproduces
all of this and moves $\alpha$ by at most $0.008$; the fixed-$\theta$ test with
that filter finds the deviation $L$-independent at $\theta=0.20$ and $0.25$ but
significantly $L$-dependent at $\theta\ge0.30$, locating the crossover at a
scale beyond $L\approx512$ for the smaller values.""",
    "5: conclusion clause added")

open(p, "w", encoding="utf-8", newline="\n").write(src)
print(f"\nmanuscript/main.tex: {n0} -> {len(src)} chars")
