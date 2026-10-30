# Task 3. FFGN mathematical audit

Date 2026-10-26. Scope is `code/check_generator_amplitude.py`, `tests/test_ffgn_covariance.py` and this note. No manuscript, simulator, figures, archived campaigns or frozen review report were edited. No new external search or full simulation campaign was performed.

## What changed

The amplitude script now evaluates the exact stationary covariance of this repository's shared-innovation AR variant. It explicitly contrasts independent innovations, retains the public functions `ffgn_seq`, `fft_ac`, `fit_amplitude`, `exact_fgn_amplitude` and `spectral_family_amplitude`, and adds the pure numerical helper `ar_sum_covariance`.

`fit_amplitude` no longer takes absolute values or selects only positive points from a signed window. Such a window returns NaN. A fully positive window still returns the geometric-mean coefficient for a fixed nominal exponent, explicitly labelled a finite-window descriptor. Main uses deterministic covariance, not empirical autocorrelation. It prints its results without replacing `data/generator_amplitude.npz` or any other old evidence file. Old plotting scripts reading that archive must not be mistaken for newly regenerated evidence.

## Exact covariance and spectrum

Let the stable finite collection of modes satisfy

\[
X_n(t)=r_nX_n(t-1)+z_n(t),\qquad
\eta(t)=\sum_{n=1}^N W_nX_n(t),\qquad |r_n|<1.
\]

Assume centred innovations independent across time, with same-time covariance
\(Q_{nm}=\operatorname{Cov}[z_n(t),z_m(t)]\). Stationarity gives

\[
\Sigma_{nm}=\operatorname{Cov}(X_n,X_m)
=r_nr_m\Sigma_{nm}+Q_{nm}
=\frac{Q_{nm}}{1-r_nr_m}.
\]

For nonnegative lag,

\[
C_N(\tau)=\operatorname{Cov}[\eta(t+\tau),\eta(t)]
=\sum_{n,m}W_nW_m r_n^\tau\frac{Q_{nm}}{1-r_nr_m}.
\]

For the current common innovation \(z_n(t)=z(t)\), every entry of Q equals \(\sigma_z^2\). All cross terms remain. For independent mode innovations, \(Q_{nm}=\sigma_z^2\delta_{nm}\), giving the distinct formula

\[
C_N^{\rm ind}(\tau)=\sigma_z^2\sum_n
\frac{W_n^2r_n^\tau}{1-r_n^2}.
\]

The equivalent stationary spectra are

\[
S_N^{\rm shared}(\omega)=\sigma_z^2
\left|\sum_n\frac{W_n}{1-r_ne^{-i\omega}}\right|^2,
\qquad
S_N^{\rm ind}(\omega)=\sigma_z^2\sum_n
\frac{W_n^2}{|1-r_ne^{-i\omega}|^2}.
\]

These equations require finite innovation variance, not Gaussianity. They describe an infinite stationary innovation history. The simulator and retained `ffgn_seq` instead initialise X to zero. With t completed innovations since that initial state, their transient covariance is

\[
C_N(t,\tau)=\sum_{n,m}W_nW_m r_n^\tau Q_{nm}
\frac{1-(r_nr_m)^t}{1-r_nr_m}.
\]

Stationary values must not silently be compared with a cold-start finite record if its slow modes have not equilibrated. No new burn-in policy or simulator change was made in this task.

Multiplying all weights by a constant q multiplies every raw covariance by q² and leaves \(C(\tau)/C(0)\) unchanged. The weights in `bd_core.ffgn_weights` omit the uniform-variance compensation factor 12 that appears in the old manuscript formula. Restoring it only as a common factor would therefore change raw variance, not the normalised covariance or sign-thresholded sequence. It cannot repair missing cross terms.

## The θ=0 limit and the limits of a power-law fit

Write the implemented weights as

\[
W_n^2=K_\theta(1-e^{-2u_n})u_n^{1-2\theta},\qquad
u_n=aB^{-n},\qquad
K_\theta=\frac{B^{1/2-\theta}-B^{\theta-1/2}}{\Gamma(2-2\theta)}.
\]

For small u, \(W(u)\simeq\sqrt{2K_\theta}\,u^{1-\theta}\). The shared-innovation process has impulse response

\[
g_j=\sum_n W_nr_n^j,\qquad
C_N(\tau)=\sigma_z^2\sum_{j\ge0}g_jg_{j+\tau}.
\]

Replacing the geometric grid by an integral gives, only in an intermediate range \(u_Nj\ll1\ll u_1j\),

\[
g_j\simeq \frac{\sqrt{2K_\theta}}{\ln B}
\int_0^\infty u^{-\theta}e^{-ju}\,du
=\frac{\sqrt{2K_\theta}\Gamma(1-\theta)}{\ln B}
j^{\theta-1}.
\]

At fixed positive \(0<\theta<1/2\), the corresponding continuum-kernel approximation yields a factor proportional to
\(B(\theta,1-2\theta)\tau^{2\theta-1}\). Its θ→0 limit is nonuniform. One must not substitute θ=0 into that large-lag coefficient.

At θ=0 the approximate kernel is c/j. The elementary identity

\[
\sum_{j=1}^{\infty}\frac{1}{j(j+\tau)}
=\frac{H_\tau}{\tau}
\simeq\frac{\ln\tau+\gamma}{\tau}
\]

then explains why a log(τ)/τ contribution is possible for the shared-innovation limit. The j=0 and short-time terms supply additional terms of order 1/τ. This is a continuum-grid mechanism argument, not an exact finite-N pure-power asymptote. A fixed geometric grid can also carry log-periodic corrections. For every fixed finite N, the exact covariance is a finite sum of decaying exponentials, so its ultimate tail is exponential. Increasing N extends an intermediate scaling range; convergence at fixed lags does not exchange the N→∞ and τ→∞ limits.

Fresh deterministic output from the revised script provides a finite-window check of this distinction. At θ=0, a=6 and B=2,

| N | 64 C(64)/C(0) | 256 C(256)/C(0) | 1024 C(1024)/C(0) | 4096 C(4096)/C(0) |
|---|---:|---:|---:|---:|
| 16 | 1.035965 | 1.231520 | 1.350184 | 1.222553 |
| 23 | 1.046737 | 1.271586 | 1.495011 | 1.714344 |
| 32 | 1.046871 | 1.272101 | 1.496967 | 1.721721 |
| 48 | 1.046871 | 1.272103 | 1.496973 | 1.721743 |

For N=48 the raw-field normalised lag-one correlations are 0.377315 for common innovations and 0.237958 for independent innovations. These are not binary-field correlations. With the automatic N=27 used for T=2²⁰, a fixed-exponent fit on τ=64–4096 gives A_window=1.563335 at θ=0, but this number must not be relabelled a constant asymptotic 1/τ amplitude. It is also not the old empirical value extracted after rectifying an autocorrelation floor.

The unbranched shared AR construction is correlated at θ=0. The production simulator has a separate θ==0 i.i.d. shortcut, so its parameter family switches construction there. This directly establishes a discontinuous input family for this implementation. It does not establish an RG phase or historical authors' exact implementation.

## Gaussianity, thresholding and physical interpretation

The implemented innovation is U[0,1]−1/2 with variance 1/12. A linear filter of uniform innovations is not strictly a Gaussian process. Approximate Gaussianity may be tested, but cannot be assumed from the label FFGN or from summing AR components that share the same innovations.

For a centred jointly Gaussian pair with normalised correlation ρ, the centred sign field s=2·1[η>0]−1 satisfies

\[
\langle s(t)s(t+\tau)\rangle=\frac2\pi\arcsin\rho(\tau).
\]

Equivalently, \(\operatorname{Cov}(b_t,b_{t+\tau})=\arcsin\rho/(2\pi)\) for b=1[η>0], with Var(b)=1/4. This identity applies to Gaussian fGn and Gaussian spectral fields. It is not exact for the uniform-driven AR variant. The latter's binary autocorrelation should be measured directly or derived from its actual joint law. Covariance of the latent non-Gaussian process alone need not determine its sign covariance.

Exact discrete fGn has H=θ+1/2 and normalised covariance at integer lag τ≥1

\[
C_{\rm fGn}(\tau)=\tfrac12[(\tau+1)^{2H}-2\tau^{2H}+(\tau-1)^{2H}],
\qquad C_{\rm fGn}(0)=1.
\]

Its asymptotic coefficient H(2H−1) vanishes at θ=0. This is a property of that particular full covariance family. The continuum KPZ problem instead specifies a low-frequency nonanalytic term and an independent short-range contribution, conventionally D0+Dθ|ω|^(−2θ). Having a different short-lag covariance or relative amplitude is not, by itself, proof of a wrong temporal exponent or a different asymptotic universality class. It can change finite-scale behaviour and crossover lengths, which remain empirical questions.

A θ_eff found by matching τ=1–4096 using 1/τ weights is a finite-window shape parameter. Changing lag range, weighting or the reference family can change it. Comparing it with 1/6 cannot assign an RG phase. In a revised manuscript it should either be removed or labelled solely as a curve-shape summary, without shading a physical short-range/long-range phase region in that coordinate.

A nearest-neighbour increment pmf in the single-step model is determined by p_flat. Its Hellinger distance may detect local changes while missing scale-dependent arrangements of steps. Fixed-θ size dependence of that local distance is a useful descriptive check, but growth of the distance is not a necessary condition for every universality crossover. Long-range scaling claims require correlations or structure factors at multiple scales, dynamic scaling, or equivalent evidence. No finite-window null result here proves a θc or that a crossover lies beyond a particular size.

## Source-grounded literature positioning

The following distinctions use only the already available first-party papers and the immutable source packet. They are not a new systematic literature search.

| Work | What was actually being identified | Boundary for the revised claim |
|---|---|---|
| Squizzato and Canet, PRE 100, 062143 (2019) | Stability of the short-range fixed point and a long-range branch in the continuum RG description. The quoted d=1 value is 1/6. | This is a continuum exponent/coupling criterion, not a finite-lag fitted-shape threshold. Local source `.lit/squizzato_canet.txt` L99–107, L1320–1360 and L1636–1705. |
| Alés and López, PRE 99, 062139 (2019) | Faceting and anomalous roughening, including a split of global and spectral roughness exponents near 0.25±0.03. | That morphological/anomalous-scaling boundary is not automatically identical to onset of any deviation from short-range KPZ. Their large sizes also differ from the present small campaign. Local source `.lit/ales_lopez2019.txt` L374–415, especially L402–405. |
| Rodríguez-Fernández et al., PRE 110, 024104 (2024) | Model- and discretisation-dependent exponents and height statistics for temporally correlated driving. | Dependence on numerical/model details was already a central subject. The present work must not claim to be the first to discover all generator or implementation sensitivity. It can contribute a precise audit of this repository's common-innovation construction and corrected production layout. Local source `.lit/fragility.txt` L601–658. |
| Jiang and Xia, CTP 78, 105602 (2026) | Height-distribution shape in discrete growth models, with a reported growth-regime marker around 0.2, compared with earlier macroscopic values around 0.25. | Shape, macroscopic exponents and saturated-distribution markers are different definitions; numerical differences do not alone prove contradictory estimates of one critical point. Local source `.lit/ae8f20.txt` L399–425 and L659–685. |
| Wang and Xia, PRE 111, 024124 (2025) | Emergence measured with Hellinger distance and novelty detection. | Those choices of diagnostics are prior art. The official abstract was checked in the review, but the full implementation, calibration and numerical threshold were not independently verified here. Do not attribute this repository's detector, run splits or fixed-significance rule to it. The immutable packet main.tex L356–372 already records this access boundary. |

The published notation ζ_t(u_n) or a named FFGN family does not establish whether all modes share an innovation stream, and does not by itself establish initialization. The available Song/Xia and Jiang/Xia formulas also need careful reading because their printed details differ from this implementation. An ambiguity in a publication is not evidence that its private code is wrong. The revised manuscript should use terms such as “the shared-innovation AR variant examined here”, with a specific formula, and reserve literature-level attribution for implementations actually verified.

FFT filtering and exact-fGn spectra are established construction methods. The supportable originality here is validation of an explicitly identified input process and its consequences in an explicitly bounded numerical experiment, together with correct uncertainty calibration. It is not invention of a spectral generator, Hellinger distance, novelty detection, or a new universal physical θ axis. No claim of first observation of a negative-θ control was independently verified in this task.

## English LaTeX passage for optional manuscript adoption

```latex
We audit the shared-innovation autoregressive construction implemented here,
without assuming that every published FFGN implementation uses the same
innovation streams. Writing
$X_n(t)=r_nX_n(t-1)+z_n(t)$ and
$Q_{nm}=\langle z_n(t)z_m(t)\rangle$, its stationary covariance is
\begin{equation}
 C_N(\tau)=\sum_{n,m}W_nW_m r_n^\tau
              \frac{Q_{nm}}{1-r_nr_m},\qquad \tau\ge0.
\end{equation}
For a common innovation, $Q_{nm}=\sigma_z^2$ and the cross terms remain;
independent innovations instead give
$C_N(\tau)=\sigma_z^2\sum_n W_n^2r_n^\tau/(1-r_n^2)$.
A common rescaling of all weights changes the variance but leaves the
normalised covariance unchanged. The simulator starts its modes at zero,
whereas this expression assumes stationarity; finite-record comparisons
must account for that distinction.

The common-innovation process has impulse response
$g_j=\sum_n W_nr_n^j$. At nominal $\theta=0$, a continuum approximation
to the geometric grid gives $g_j\propto j^{-1}$ over an intermediate
range. The identity
$\sum_{j=1}^{\infty}[j(j+\tau)]^{-1}=H_\tau/\tau$
then permits a logarithmic correction to a $1/\tau$ tail. This argument
does not establish a pure power law at finite $N$: the exact finite-mode
covariance ultimately decays exponentially. We therefore report exact
finite-$N$ correlations and label fixed-exponent fits only as
finite-window summaries.

The innovations used here are centred uniform variables, so the resulting
field is not strictly Gaussian. The relation
$C_{\mathrm{sign}}=(2/\pi)\arcsin\rho$ is exact for a centred jointly
Gaussian field and is not imposed on this non-Gaussian construction.
Likewise, a parameter obtained by fitting a finite-lag correlation curve
is not identified with the temporal exponent of the continuum RG theory.
Our comparisons establish sensitivity to the specified input processes
within the simulated window, rather than a recalibration of an RG phase
boundary or a determination of a published crossover threshold.
```

## Verification record

1. RED was observed before implementation using

   `& 'D:/TOOLS/anaconda3/python.exe' -m pytest tests/test_ffgn_covariance.py -q`

   Result was 8 failed, 1 passed, exit 1. Seven failures were the intentionally absent covariance helper; the signed-floor regression failed on the old absolute-value fitter, which returned 0.024042391576364685 for an alternating ±0.01 floor. The positive-power-law regression already passed.

2. GREEN after implementation used the same command. All 9 parametrised cases passed, exit 0. They cover one AR mode in both innovation configurations, two-mode closed forms with and without cross terms, a direct 100-term impulse-response sum, raw q² scaling and normalised invariance under a common negative weight scale, and the signed-floor/positive-power-law fits.

3. `& 'D:/TOOLS/anaconda3/python.exe' code/check_generator_amplitude.py` completed with exit 0. The deterministic values quoted above came from this run. It did not write old or new amplitude archives.

4. A full `python -m pytest tests -q` run initially reached 33 passed and 3 setup errors. The affected tests were `test_spectral_simulation_archive_keeps_layout_version[powerlaw]`, the `[fgn]` instance, and `test_legacy_archive_is_not_labelled_as_corrected`. All three errors were an OSError creating pytest's system-Temp directory, before their assertions ran.

5. Retried the entire currently available suite using an explicit writable workspace temporary directory

   `& 'D:/TOOLS/anaconda3/python.exe' -m pytest tests -q --basetemp 'E:/应用物理学/科研/BD模型/docs/revision_20261026/task3_pytest_tmp' --tb=short`

   Result was 36 passed, exit 0. The temporary directory contains only pytest test artifacts. This is a snapshot of the concurrently edited project, not a substitute for the root agent's final integration test after all tasks finish.

No git commit was made. The frozen reviewer_3.md was not changed.
