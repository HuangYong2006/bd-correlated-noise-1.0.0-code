# Numerical side of KPZ / discrete growth with power-law temporal noise correlations (θ)

**Scope.** Extraction of verified, quantitatively specific facts from the numerical literature on surface growth
(KPZ equation and discrete models: ballistic deposition BD, restricted solid-on-solid RSOS) driven by noise whose
temporal correlations are power-law, `⟨η(x,t)η(x',t')⟩ = 2D δ(x-x') |t-t'|^{2θ-1}` (equivalently `D(ω) ~ ω^{-2θ}`),
`θ ∈ (0, 1/2)`. Theory threshold under test: `θ_th = 1/6` in d=1 (Medina et al. 1989; Squizzato & Canet 2019).

**Method.** All numeric claims below were read out of full texts I downloaded and extracted locally (PyMuPDF), or
from Crossref / OpenAlex / PubMed E-utilities metadata. Every claim carries a URL. Where I could not get the full
text I say so explicitly and do **not** infer numbers.

Local extracted texts (workspace, for re-checking): `.lit\2102.13356.txt`, `.lit\ae8f20.txt`, `.lit\fragility.txt`,
`.lit\ales_lopez2019.txt`, `.lit\squizzato_canet.txt`, `.lit\hu_hao_xia.txt`, `.lit\ctp_adfa78.txt`.

---

## 0. CORRECTIONS TO THE BRIEF (these matter for citation accuracy)

| Brief said | Verified record | Evidence |
|---|---|---|
| Song & Xia, J. Stat. Mech. (2021) **063213** | Song & Xia, *J. Stat. Mech.* **2021**, **073203** (DOI 10.1088/1742-5468/ac06c3) | Crossref [api.crossref.org/works/10.1088/1742-5468/ac06c3](https://api.crossref.org/works/10.1088/1742-5468/ac06c3) |
| "Squizzato, Canet et al., *Extensive numerical simulations of surface growth with temporally correlated noise*, PRE **103**, 012121 (2021)" | **Mis-attributed.** PRE 103, 012121 = **Song, T. & Xia, H.**, *Extensive numerical simulations of surface growth with temporally correlated noise*. It is a **ballistic-deposition** study, not Squizzato/Canet. | Crossref [api.crossref.org/works/10.1103/PhysRevE.103.012121](https://api.crossref.org/works/10.1103/PhysRevE.103.012121); PubMed [33601554](https://pubmed.ncbi.nlm.nih.gov/33601554/) |
| Squizzato & Canet 2019 = numerical | Squizzato & Canet, *KPZ equation with temporally correlated noise: A nonperturbative renormalization group approach*, PRE **100**, 062143 (2019), arXiv:1907.02256. **Analytic (NPFRG) only — no numerical simulations.** | Crossref [10.1103/PhysRevE.100.062143](https://api.crossref.org/works/10.1103/PhysRevE.100.062143); [arxiv.org/abs/1907.02256](https://arxiv.org/abs/1907.02256) |
| "Height distributions of discrete growth models…, Commun. Theor. Phys. (2025)" | Published as *Commun. Theor. Phys.* **78**, 105602 (**2026**), received 5 Feb 2026, published 7 Sep 2026 (per the article PDF/publisher record). Open access, CC BY-NC-ND. | [iopscience.iop.org/article/10.1088/1572-9494/ae8f20](https://iopscience.iop.org/article/10.1088/1572-9494/ae8f20) |

**Theory anchors (verified wording/values):**
- Medina, Hwa, Kardar & Zhang 1989 (PRA 39, 3053), as reproduced by Rodríguez-Fernández et al. 2024:
  `α = 1/2 (θ<1/6)`, `α = 1.69θ + 0.22 (1/6<θ<1/2)`; `z = 3/2 (θ<1/6)`, `z = (2α+1)/(1+2θ) (θ>1/6)`;
  "no stable fixed point is found for θ > 1/4… an infinite set of singularities appear at θ* = 1/4, 2/6, 3/8, …".
  Source: [PRE 110, 024104 (2024)](https://repositorio.unican.es/xmlui/bitstream/10902/33783/2/FragilityKardarParisi.pdf), Eqs. (3)–(4) and text.
- Squizzato & Canet 2019: **"the SR-KPZ fixed point is stable below a threshold θ_th = 1/6 in d = 1 and θ_th ≃ 0.35
  in d = 2"**; their flow gives "θ_th = 0.166 confirming the theoretical prediction θ_th = 1/6".
  Source: [arXiv:1907.02256](https://arxiv.org/abs/1907.02256) full text (`squizzato_canet.txt`).
- Competing analytic results with **no threshold**: Katzav & Schwartz SCE (PRE 70, 011601 (2004)) and Fedorenko FRG
  (PRB 77, 094203 (2008)); "for any non-zero θ (no threshold value)".
  Source: [arXiv:1907.02256](https://arxiv.org/abs/1907.02256); [PRE 110, 024104 (2024)](https://repositorio.unican.es/xmlui/bitstream/10902/33783/2/FragilityKardarParisi.pdf).

---

## 1. Song & Xia — "Kinetic roughening and nontrivial scaling in the KPZ growth with long-range temporal correlations"

*J. Stat. Mech.* **2021**, 073203 (DOI [10.1088/1742-5468/ac06c3](https://doi.org/10.1088/1742-5468/ac06c3));
arXiv:[2102.13356v2](https://arxiv.org/abs/2102.13356v2). **Full text obtained** (arXiv v2 PDF, 30 pp).

### 1a. Exact noise-generation recipe and run parameters
- **Method: fast fractional Gaussian noise (FFGN), Mandelbrot generator.** Not Fourier filtering, and **no**
  rank/quantile mapping. They state explicitly that the alternative Fourier filtering method "requires generating
  the complete temporal sequence of random numbers for each site, which means FFM takes up more memory and time
  than the FFGN" (comparison is in the companion paper [arXiv:2105.13528](https://arxiv.org/abs/2105.13528)).
- The imposed quantity is the **time autocorrelation**, not a spectral exponent directly:
  `C(τ) ≡ ⟨η(t)η(t+τ)⟩ ~ τ^{2θ−1} (τ→∞)` — equivalent to `D(ω) ~ ω^{−2θ}`.
- FFGN construction (their Eq. 4): `η(t) = Σ_{n=1}^{N} W_n X_t(u_n)`, with `u_n = a B^{−n}`, `r_n = e^{−u_n}`,
  `B = 2`, `a = 6`, `ζ(t)` uniform on [0,1],
  `W_n² = [12(1−r_n²)/Γ(2−2θ)] (B^{1/2−θ} − B^{θ−1/2}) u_n^{1−2θ}`,
  `X_1 = (1−r_n²)^{−1/2}(ζ(1)−1/2)`, `X_t = r_n X_{t−1} + (ζ(t)−1/2)`.
- **Number of components: `N` chosen in the range [20, 45]** "relying on the sequence length we used".
- `{η(t)}` is an **independent Mandelbrot sequence for every fixed lattice site** (spatially uncorrelated).
- **No binarisation / no integer mapping** in this paper. They explicitly attribute the discrepancy between their
  continuum KPZ results and the discrete BD results to the fact that "one needs binarizing the correlated noise
  into integers in the modified BD model, while [this] does not introduce this binarization in the discretized
  version of the temporal correlated KPZ equation."
- **Integration scheme:** nonlinear term replaced by an exponentially decaying function
  `f(x) = (1 − e^{−cx})/c` (Dasgupta et al. stabilisation), `ν = 1`, `λ = 4`, `c = 0.1`, `Δx = 1`, `Δt = 0.05`,
  flat initial condition `h(x,0)=0`, periodic boundary conditions. Hybrid GPU+CPU runs.
- **(1+1) d:** `L = 65536` for `W(L,t)` (Fig. 2); `L = 32768` for `G(l,t)` and `S(k,t)` (Figs. 4–5) with
  local window `2 ≤ l ≤ 1024` and structure-factor times `t = 1.28×10³, 5.12×10³, 2.05×10⁴, 8.19×10⁴, 3.28×10⁵`;
  morphology snapshots at early `t = 1.28×10³` and saturated `t = 4.81×10⁶`. Averaging: 500 noise realizations
  (Figs. 2, 4); 200 (θ=0.20) and 150 (θ=0.40) realizations (Fig. 5).
- **(2+1) d:** `L × L = 4096 × 4096` for `G(l,t)` with `2 ≤ l ≤ 128`, at `t = 4.10×10⁴`; and `2048 × 2048` for
  `S(k,t)` with times `1.60×10², 6.40×10², 2.56×10³, 1.02×10⁴, 4.10×10⁴`; 150 realizations; morphology at early
  `t = 6.3×10²` and saturated `t = 1.30×10⁶`.

### 1b. Measured exponents (only these two θ values are given numerically in the text)
**(1+1) d**
| θ | β | α (from W(L,t)) | α (from S(k,t) collapse) | α_s | z |
|---|---|---|---|---|---|
| 0.20 | **0.41** | **0.60** | 0.61 | **0.66** | 1.47 (1.45 from S) |
| 0.40 | **0.65** | **0.90** | 0.90 | **0.95** | 1.37 |

**(2+1) d**
| θ | β | α (from W) | α (from S collapse) | α_s | z |
|---|---|---|---|---|---|
| 0.20 | **0.28** | **0.48** | 0.50 | **0.54** | 1.73 (1.68 from S) |
| 0.40 | **0.59** | **0.82** | 0.88 | **0.96** | 1.38 (1.36 from S) |

The full θ-dependence of α, β, z is presented only graphically (their Fig. 6 for 1+1 d, Fig. 12 for 2+1 d), as
comparisons against DRG (Medina), SCE (Katzav–Schwartz), FRG (Fedorenko) and NPFRG (Squizzato–Canet) curves.
I could not read numeric values off those figures.

### 1c. Threshold conclusion — explicit
- They **do** claim a critical threshold exists: "there should exist a critical threshold θ_th separating the small
  and large θ regimes, as one can observe from θ = 0.20 to θ = 0.40 in Figure 7(f)-(g)"; no obvious mounds at
  θ = 0.05 and 0.20, obvious mounds at θ = 0.40, full faceting at θ = 0.48.
- They **do not** locate it at 1/6. They cite the value **θ_th ≈ 1/4** from Refs. [11,33] (= Alés & López 2019,
  PRE 99, 062139; Alés & López 2020, J. Stat. Mech. 033210) and conclude: "the critical threshold in the stochastic
  growth model driven by the temporally correlated noise is **universal, which does not rely on the substrate
  dimension**."
- Crucially they also assert **nontrivial scaling at small θ**, which is the opposite of the Medina/Squizzato–Canet
  picture: "our results provide numerical evidences that **nontrivial scaling can occur even in the small θ
  regime**… which is a little different with the recent results [11]" (i.e. different from Alés & López's
  θ_th ≈ 0.25 standard-scaling window).
- For θ = 0.20 > 1/6 they get α = 0.60, β = 0.41 — **not** the SR-KPZ values (α=1/2, β=1/3) and not Medina's
  one-loop `α = 1.69(0.2)+0.22 = 0.558` either.

### 1d. Faceted morphology and how it is quantified
- "when θ is close to 0.5, the temporal correlated KPZ interfaces develop full faceted patterns in the long growth
  time limits"; "mounded structures on different scales"; "the obvious segmentation happens when the temporal
  correlation exponent is beyond a certain critical threshold".
- **Quantification is indirect, not a facet-detection metric:** (i) by-eye interface profiles, (ii) the observation
  that "when θ is near 0.5, the critical values of both α and α_s approach to 1", (iii) the appearance of
  `α_s > α` splitting. No facet-size distribution, no order parameter, no clustering metric is used.

### 1e. Anomalous scaling — explicit
- **Yes.** Abstract: "the characteristic roughness exponents satisfy `α < α_s`", and "α_loc exhibiting
  non-universal scaling within local window sizes". Text: "We find that the global roughness exponents are less
  than the corresponding spectral scaling exponents… Thus our results provide numerical evidences that nontrivial
  scaling can occur even in the small θ regime" and "It is remarkable that `α < α_s` does not satisfy the existing
  dynamic scaling classifications."
- **α_loc is not a universal number:** "the height–height correlation function does not exhibit universal scaling
  behavior within local window sizes, and thus the local roughness exponent is not universal, more precisely, the
  local roughness exponents tend to decrease with window size. Therefore, there are no evidences for the existence
  of a universal critical exponent describing local window". Same statement is made for (2+1) d.
- So: **α ≠ α_s (global vs spectral)** is claimed, and **α_loc is explicitly non-universal / window-dependent**
  rather than a well-defined second exponent.

---

## 2. Song & Xia, Phys. Rev. E 103, 012121 (2021) — "Extensive numerical simulations of surface growth with temporally correlated noise"

DOI [10.1103/PhysRevE.103.012121](https://doi.org/10.1103/PhysRevE.103.012121); PubMed
[33601554](https://pubmed.ncbi.nlm.nih.gov/33601554/). **PAPER IS PAYWALLED (APS, `oa_status: closed` per Unpaywall);
no arXiv version exists** (arXiv API title queries return nothing; the arXiv `cond-mat.stat-mech` listing of papers
with "temporally correlated noise" contains no such entry). **I could not access the full text and therefore cannot
report β(θ) or α(θ) numbers.**

**Verified abstract (verbatim, PubMed E-utilities):**
> "Surface growth processes can be significantly affected by long-range temporal correlations. In this work, we
> perform extensive numerical simulations of a (1+1)- and (2+1)-dimensional **ballistic deposition (BD) model**
> driven by temporally correlated noise, which is regarded as the temporal correlated Kardar-Parisi-Zhang
> universality class. Our results are compared with the existing theoretical predictions and numerical simulations.
> **When the temporal correlation exponent is above a certain threshold, BD surfaces develop gradually faceted
> patterns.** We find that the temporal correlated BD system displays nontrivial dynamic properties, and the
> characteristic roughness exponents satisfy **α ≃ α_loc < α_s in (1+1) dimensions**, which is beyond the existing
> dynamic scaling classifications."

**Model definition (verified from secondary full texts, not from the paper itself):**
- Jiang & Xia (2026), [doi:10.1088/1572-9494/ae8f20](https://iopscience.iop.org/article/10.1088/1572-9494/ae8f20),
  Eq. (10), attribute the "temporally correlated BD model [48, 60, 62]" rule to refs [48] = this paper,
  [60] = Wang & Xia PRE 111 024124, [62] = Lam, Sander & Wolf PRA 46 R6128 — namely
  `h(x,t+1) = max[h(x−1,t), h(x,t) + η'_θ(x,t), h(x+1,t)]` with binarised noise `η' = 1 if η>0 else 0`.
- Rodríguez-Fernández et al. (2024) describe Song & Xia (2021a) as "simulations of BD with correlated noise,
  generated by two different algorithms (one of them was identical to that of Lam et al.)" and that they
  "concluded that the prediction of the DRG was in excellent agreement with simulations"
  ([Fragility, PRE 110, 024104](https://repositorio.unican.es/xmlui/bitstream/10902/33783/2/FragilityKardarParisi.pdf)).
- Song & Xia (2021b, JSTAT 073203) note their own continuum values are "somewhat **larger** than the previous
  results of the BD model in presence of long-range temporal correlations [13,14]" where [13] = this paper.
  No numbers are quoted.

**Threshold:** the abstract states a threshold exists ("above a certain threshold… faceted patterns") but gives no
value. **I could not verify any numeric θ_th from this paper.**

Note: `α ≃ α_loc < α_s` here is *different* from the continuum JSTAT 2021 result (`α_loc` non-universal, `α < α_s`),
i.e. the discrete BD and continuum KPZ numerics of the same group differ on the status of `α_loc`.

---

## 3. Song & Xia, J. Stat. Mech. (2016) 113206 — "Long-range temporal correlations in the KPZ growth: numerical simulations"

DOI [10.1088/1742-5468/2016/11/113206](https://doi.org/10.1088/1742-5468/2016/11/113206). **Paywalled
(`oa_status: closed`); no arXiv version found; Semantic Scholar reports the abstract as elided by the publisher.
I could not access the full text, so I cannot report the model definition in detail nor any measured exponents.**

**Verified abstract (OpenAlex, [api.openalex.org/works/doi:10.1088/1742-5468/2016/11/113206](https://api.openalex.org/works/doi:10.1088/1742-5468/2016/11/113206)):**
> "To analyze long-range temporal correlations in surface growth, we study numerically the **(1+1)-dimensional
> Kardar–Parisi–Zhang (KPZ) equation driven by temporally correlated noise**, and obtain the scaling exponents
> based on **two different numerical methods**. Our simulations show that **the numerical results are in good
> agreement with the dynamic renormalization group (DRG) predictions**, and are also consistent with the
> simulation results of the ballistic deposition (BD) model."

**Conflicting secondary characterisation (report as-is, unresolved):** Rodríguez-Fernández et al. 2024 describe
this paper as "simulations of **ballistic deposition** with correlated noise, generated by two different algorithms
(one of them was identical to that of Lam et al. [25])"
([Fragility](https://repositorio.unican.es/xmlui/bitstream/10902/33783/2/FragilityKardarParisi.pdf)), whereas the
abstract says the **KPZ equation**. I cannot resolve which is right without the full text.

---

## 4. Wang & Xia, "Emergence in kinetic roughening with long-range temporal correlations", Phys. Rev. E 111, 024124 (2025)

DOI [10.1103/PhysRevE.111.024124](https://doi.org/10.1103/PhysRevE.111.024124); PubMed
[40103081](https://pubmed.ncbi.nlm.nih.gov/40103081/); Crossref authors **Wang, Shuting; Xia, Hui**.
**Paywalled (Unpaywall `oa_status: closed`; no arXiv preprint found). Full text NOT accessible.**

**Verified abstract (verbatim, PubMed E-utilities / OpenAlex):**
> "The role of long-range temporal correlations in kinetic roughening processes is significant, potentially
> influencing surface morphologies and dynamic scaling properties. This type of correlation, acting as a nonlocal
> interaction during growth processes, disrupts the self-affine fractal structure by causing the emergence of
> nontrivial global properties from individual local interactions. In this work, **two typical discrete growth
> models, including random deposition and ballistic deposition, are investigated numerically in the presence of
> long-range temporal correlations. Two quantitative measurement methods, the Hellinger distance and the novelty
> detection, are used independently to determine whether emergence occurs.** Our simulation results indicate that,
> **above a critical threshold of temporal correlation exponent, the systems exhibit emergence behavior**,
> consistent with the existence of nontrivial scaling properties and surface morphologies when long-range temporal
> correlations are presented."

**Model:** the two models are **random deposition (RD)** and **ballistic deposition (BD)** — *not* RSOS and not the
continuum KPZ equation.

**Threshold value:** not in the abstract. Jiang & Xia (2026) state it explicitly:
"studies focusing on dynamic properties, including anomalous scaling [46, 59] and **emergent phenomena [60], have
reported a critical value θ = 0.25**" ([doi:10.1088/1572-9494/ae8f20](https://iopscience.iop.org/article/10.1088/1572-9494/ae8f20),
where [60] = this paper). This is **θ_c = 0.25**, i.e. **not** 1/6.

**Does it argue for "emergence" of a new universality class?** It frames the finding as an **emergence
phenomenon** (nontrivial global behaviour from local interactions) and states that the threshold behaviour is
"consistent with the existence of nontrivial scaling properties and surface morphologies". Whether the wording
"new universality class" appears in the body **could not be verified** (full text inaccessible).

**What I could NOT verify (important):** the exact novelty-detection algorithm (Isolation Forest? One-Class SVM?
Local Outlier Factor? Autoencoder reconstruction error?), the features used, the reference ensemble for the
Hellinger distance, and the system sizes/times. These are not in the abstract and the full text is paywalled.

---

## 5. Jiang & Xia, "Height distributions of discrete growth models in the presence of long-range temporal correlations"

*Commun. Theor. Phys.* **78**, 105602 (2026); DOI [10.1088/1572-9494/ae8f20](https://doi.org/10.1088/1572-9494/ae8f20).
**Open access (CC BY-NC-ND). Full text PDF obtained from
`https://beta.iopscience.iop.org/article/10.1088/1572-9494/ae8f20/pdf` and extracted.**

### 5a. Exact model definition and noise recipe
- Noise autocorrelation imposed: `⟨η(x,t)η(x′,t′)⟩ ~ δ^d(x−x′)|t−t′|^{2θ−1}`, `θ ∈ [0, 0.4]`. They restrict to
  θ ≤ 0.4 because "the FFGN method becomes numerically unstable for large temporal correlation exponents θ > 0.4".
- **FFGN (Mandelbrot) + binary mapping.** Weights (their Eq. 6):
  `W_n² = [12/Γ(2−2θ)](1−r_n²)(B^{1/2−θ} − B^{θ−1/2})(aB^{−n})^{1−2θ}` with `B = 2`, `a = 6`, `u_n = aB^{−n}`,
  `r_n = e^{−u_n}`. Recursion (Eq. 7): `X_1 = |ζ_1−1/2|/√(1−r²)`, `X_t = r X_{t−1} + |ζ_t−1/2|`.
  Sum (Eq. 8): `η_θ(t) = Σ_{n=1}^{N} W_n X_t(u_n)` with **N = 40**.
- **Binary mapping (Eq. 9): `η'_θ(t) = 1 if η_θ(t) > 0, else 0`** (threshold at zero; attributed to
  Lam, Sander & Wolf, PRA 46 R6128 (1992)).
- **BD rule (Eq. 3, standard):** `h(x,t+1) = max[h(x−1,t), h(x,t) + 1, h(x+1,t)]` — the **true ballistic rule**
  (allows |Δh| > 1).
- **Correlated BD (Eq. 10):** `h(x,t+1) = max[h(x−1,t), h(x,t) + η'_θ(x,t), h(x+1,t)]`.
- **RSOS (Eqs. 4, 11):** `h(x,t+1) = h(x,t) + η'_θ(x,t)` if `|h(x,t) − h(x±1,t)| ≤ 1`, else `h(x,t+1)=h(x,t)`.
- Flat initial condition, periodic boundaries (leftmost/rightmost columns nearest neighbours), parallel update.
- **Sizes/times:** substrate size **L = 1024** ("to suppress finite-size effects"); early-time/intermediate/saturation
  snapshots at `t = 200`, `t = 2000`, and saturation `t = 218` as printed (this last value is almost certainly a
  typesetting/PDF-extraction artefact of a power of ten or of `2^18`; I report it verbatim rather than guess).
  Characteristic analysis times are `t_1` (growth regime) and `t_sat`, both read off the `W(t)` curves, not given
  as numbers.
- RSOS is compared after the inversion `h → −h`, because BD ↔ λ>0 and RSOS ↔ λ<0.

### 5b. How the height distribution is analysed — the critical question
- **They compare the rescaled distributions directly to Tracy–Widom distributions.** Figures 2 (BD) and 4 (RSOS)
  are log-linear plots of the height distribution, rescaled to zero mean and unit variance, "compared with three
  precisely rescaled TW distributions" (i.e. **GOE, GUE, GSE**), with insets zooming on the peaks. Saturation-regime
  height-increment distributions are compared to the **Baik–Rains BR-F0** distribution and to a Gaussian.
- **No distribution-distance metric is used.** I grepped the full extracted text for `Hellinger`, `Kolmogorov`,
  `Kullback`, `total variation`, `fidelity`, `chi-square`, `Isolation`, `forest`, `machine learning`, `novelty`:
  **zero occurrences of every one of these terms.**
- The **only quantitative** measures are the **skewness S (Eq. 12) and kurtosis K (Eq. 13)**, compared against the
  **exact TW-GOE values quoted in the paper: `S = 0.2934`, `K = 3.1652`** (Fig. 3 and Fig. 5 solid lines) and the
  **BR-F0 values `S = 0.3594`, `K = 3.2891`** (Fig. 9 and Fig. 11 solid lines). Gaussian references `S = 0, K = 3`.
  S and K are plotted vs θ with error bars.
- So the "distance to the reference ensemble" is a **two-moment comparison**, not a full distribution distance.

### 5c. Conclusions about whether the height distribution stays KPZ/TW-GOE
- **Growth regime, BD:** "For BD within a small θ region, namely 0 ≤ θ ≤ 0.2, the effects of long-range temporal
  correlations on the height distribution are very weak, and all height distributions are qualitatively similar to
  the TW-GOE form… for θ ≥ 0.25, the right tail exhibits a progressive deviation, and the peak value gradually
  shifts upward, indicating a **continuous crossover from the KPZ universality class to an unknown universality
  class**." Estimated S and K "remain stable and approach the exact values for TW-GOE up to θ = 0.2"; for θ > 0.2
  both "exhibit a pronounced increase".
- **Growth regime, RSOS:** same qualitative behaviour; "for small values of θ… distributions are closely similar to
  the TW-GOE distribution form"; pronounced increase of S,K in `0.2 ≤ θ ≤ 0.4`.
- **Critical threshold: `θ_c ≈ 0.2`** for both BD and RSOS from height distributions. Explicit: "We identify a
  critical threshold, θ_c ≈ 0.2, which separates two distinct regimes: for θ < θ_c, the system exhibits a stable
  plateau consistent with the TW-GOE distribution, indicating the robustness of the KPZ universality class.
  Conversely, for θ > θ_c, the observed deviation signifies a crossover toward a new universality class."
- They note this differs from the macroscopic-dynamics threshold: "this value differs from the threshold θ = 0.25
  derived from macroscopic dynamic phenomena [46, 60], suggesting that height fluctuations are more sensitive to
  long-range temporal correlations."
- **Saturation regime:** "the height-increment distributions of both the standard BD and RSOS models converge to
  the BR-F0 distribution." Under temporal correlations the BD increment distribution's right tail shifts left and
  the peak is suppressed for θ ∈ [0, 0.3], then for 0.3 < θ < 0.4 it becomes right-skewed with a rising peak;
  RSOS shows the same non-monotonic behaviour with the crossover at **θ = 0.25**.
- Overall: "the two models with long-range temporal correlations belong to a **continuously varying universality
  class** that depends on the temporal correlation exponent", and "a precise theoretical prediction for the KPZ
  universality class with long-range temporal correlations is still lacking."

---

## 6. Machine learning in kinetic roughening / surface-growth universality classification

### 6a. Verified affirmative findings
1. **Wang & Xia, PRE 111, 024124 (2025)** — see §4. It is **published prior art for applying a
   distribution-distance metric (Hellinger distance) and "novelty detection" to locate a threshold/crossover in
   exactly this model class (RD and BD with power-law temporal noise)**. Full text paywalled; **the specific
   novelty-detection algorithm could not be verified.** Nothing in the accessible record says "Isolation Forest".
2. **Public GitHub repository** `adamfbentley/ml-universality-classification` —
   "Exploratory anomaly detection for simulated surface-growth universality classes"
   ([github.com/adamfbentley/ml-universality-classification](https://github.com/adamfbentley/ml-universality-classification);
   README via [raw.githubusercontent.com/.../README.md](https://raw.githubusercontent.com/adamfbentley/ml-universality-classification/main/README.md)).
   Verified README content: trains on simulated EW and KPZ trajectories; includes "**Isolation Forest, Local
   Outlier Factor, and One-Class SVM comparisons**"; 16-feature pipeline; test dynamics include MBE, VLDS, quenched
   KPZ and ballistic deposition; explicitly self-described as "an exploratory computational study, **not** a claim
   of a new universal invariant" and "Research prototype". **This is not peer-reviewed literature.** I could not
   verify authorship or its relationship to any published paper.

### 6b. Verified negative findings (arXiv full-text field search, 2026-09-21)
Using the arXiv API (`export.arxiv.org/api/query`), `all:` field searches:
- `"machine learning" AND "kinetic roughening"` → **0 results**
  ([query](https://export.arxiv.org/api/query?search_query=all:%22machine+learning%22+AND+all:%22kinetic+roughening%22)).
- `"ballistic deposition" AND "machine learning"` → **0 results**
  ([query](https://export.arxiv.org/api/query?search_query=all:%22ballistic+deposition%22+AND+all:%22machine+learning%22)).
- `"surface growth" AND "anomaly detection"` → **0 results**
  ([query](https://export.arxiv.org/api/query?search_query=all:%22surface+growth%22+AND+all:%22anomaly+detection%22)).
- `"isolation forest" AND "universality"` → 5 results, **all computer-science/engineering** (LLM agent defences,
  laughter segmentation, jailbreak detection, geared-motor acoustics, bearing fault detection) — **none in physics
  or statistical mechanics** ([query](https://export.arxiv.org/api/query?search_query=all:%22isolation+forest%22+AND+all:%22universality%22)).
- General web searches for "machine learning universality class surface growth", "unsupervised learning crossover
  kinetic roughening", "Hellinger distance universality class growth", "isolation forest phase transition
  statistical mechanics" surfaced **no peer-reviewed paper** applying ML/anomaly detection to kinetic-roughening
  universality classification.

---

## 7. Which discrete update rule does each paper actually use?

**Canonical rule and its variants (verified quotations):**
- **True ballistic deposition (BD)** — `h_i → max(h_{i−1}, h_i + 1, h_{i+1})`, deposits at the first-contact
  nearest-neighbour occupied site, allows height differences > 1:
  Jiang & Xia 2026, Eq. (3): `h(x,t+1) = max[h(x−1,t), h(x,t) + 1, h(x+1,t)]`, with the physical description
  "all incident particles are released from a randomly selected position far above the aggregate, follow a vertical
  trajectory until they reach the aggregate front, and stick upon the first contacting nearest neighbor occupied
  site". [doi:10.1088/1572-9494/ae8f20](https://iopscience.iop.org/article/10.1088/1572-9494/ae8f20)
- **RSOS** — `h_i → h_i + 1` only if all `|Δh| ≤ 1`, otherwise the attempt is rejected (equivalent to
  `h_i → max(h_{i−1}, h_{i+1}) + 1` implemented via rejection): Jiang & Xia 2026, Eq. (4). Same URL.

| Paper | Discrete model | Update rule actually used | Evidence |
|---|---|---|---|
| **Jiang & Xia 2026** (CTP 78, 105602) | standard BD | `max[h(x−1,t), h(x,t)+1, h(x+1,t)]` — **true ballistic** | Eq. (3), full text |
| **Jiang & Xia 2026** | **BD with temporally correlated deposition** | `max[h(x−1,t), h(x,t)+η'_θ(x,t), h(x+1,t)]` — **true ballistic, binarised noise** (not RSOS) | Eq. (10), full text |
| **Jiang & Xia 2026** | RSOS with temporally correlated noise | `h(x,t)+η'_θ` subject to `|Δh| ≤ 1`, else rejected — **RSOS rule, kept as a separate model** | Eq. (11), full text |
| **Alés & López 2019** (PRE 99, 062139; arXiv:1902.06674) | BD with temporally correlated noise | `h(x,t+1) = Max[h(x,t) + η(x,t), h(x−1,t), h(x+1,t)]`, η ∈ {0,1} — **true ballistic**; "This digitalization of the noise enhances the statistics … by avoiding the formation of overhangs on the surface" | [arXiv:1902.06674](https://arxiv.org/abs/1902.06674) full text |
| **Rodríguez-Fernández et al. 2024** (PRE 110, 024104) "BD-I" | BD, Lam et al. noise | `h_j(t+1) = Max[h_j(t) + ζ_j(t), h_{j−1}(t), h_{j+1}(t)]`, ζ∈{0,1}, `ζ=1 iff η>0` — **true ballistic** | [Fragility full text](https://repositorio.unican.es/xmlui/bitstream/10902/33783/2/FragilityKardarParisi.pdf), §II C |
| **Rodríguez-Fernández et al. 2024** "BD-II" | BD, integer noise | same rule with `ζ = [[η]] + 1` for η > 0, else 0 — **true ballistic** | same |
| **Wang & Xia 2025** (PRE 111, 024124) | RD + BD | not verifiable (paywalled); Jiang & Xia 2026 ref. [60] groups it with the `max[h(x−1,t), h(x,t)+η', h(x+1,t)]` correlated-BD rule | attribution in CTP Eq. (10) |
| **Song & Xia 2021a** (PRE 103, 012121) | BD | not verifiable directly (paywalled); Jiang & Xia 2026 ref. [48] attributes the correlated-BD `max[…, h+η', …]` rule to it | attribution in CTP Eq. (10) |
| **Song & Xia 2016** (JSTAT 113206) | KPZ equation per abstract (BD per Fragility) | not verifiable (paywalled) | abstract text in §3 |
| **Song & Xia 2021b** (JSTAT 073203) | continuum KPZ equation | no discrete rule; noise **not** binarised | full text §2 |

**Answer to the question posed:** the discrete "BD with temporally correlated deposition" in this literature uses
the **true ballistic rule** `h_i → max(h_{i−1}, h_i + η'_i, h_{i+1})` with a **binarised (0/1) temporally
correlated noise**, i.e. height differences > 1 are allowed. The **RSOS rule** is used only where RSOS is
explicitly the model (Jiang & Xia 2026), not as the BD update. The binarisation is itself recognised in the
literature as a source of discrepancy between discrete BD and continuum KPZ results (Song & Xia 2021b,
Conclusions: "one needs binarizing the correlated noise into integers in the modified BD model, while [this] does
not introduce this binarization in the discretized version of the temporal correlated KPZ equation").

---

## 8. Consolidated table: θ vs measured (β, α, α_s, z)

All entries are **numerical measurements** unless the row is marked theory. "—" = not reported / not accessible.

| System (paper) | θ | β | α | α_s | z | Notes |
|---|---|---|---|---|---|---|
| continuum KPZ (1+1) d — Song & Xia, JSTAT 2021 (073203) | 0.20 | 0.41 | 0.60 | 0.66 | 1.47 | L=65536; 500 realizations |
| continuum KPZ (1+1) d — Song & Xia, JSTAT 2021 | 0.40 | 0.65 | 0.90 | 0.95 | 1.37 | same |
| continuum KPZ (2+1) d — Song & Xia, JSTAT 2021 | 0.20 | 0.28 | 0.48 | 0.54 | 1.73 | 150 realizations |
| continuum KPZ (2+1) d — Song & Xia, JSTAT 2021 | 0.40 | 0.59 | 0.82 | 0.96 | 1.38 | same |
| continuum KPZ (1+1) d — Hu, Hao & Xia, Physica A 619, 128744 (2023) | 0.30 | 0.583 ± 0.008 (PS) / 0.580 ± 0.010 (LS) | 0.73 ± 0.03 (both) | 0.828 ± 0.006 (PS) / 0.827 ± 0.006 (LS) | — | L=2048…16384; 750 realizations; extrapolations α_s(t) ≈ α_s + C t^{−Δ3}, α_loc with Δ2 = 0.6; DRG predicts α ≈ 0.727 at θ=0.30; "±" is rendered "+" in the PDF text layer |
| **BD (1+1) d — Alés & López, PRE 99, 062139 (2019)** | 0.15 | — | 0.56 ± 0.02 | 0.56 ± 0.02 | 1.50 | L = 2048/4096/8192; 100 realizations |
| **BD (1+1) d — Alés & López 2019** | 0.47 | — | 1.02 ± 0.01 | 1.23 ± 0.02 | 1.50 | faceted phase |
| **BD (1+1) d — Song & Xia, PRE 103, 012121 (2021)** | any θ | **not accessible** | **not accessible** | **not accessible** | — | abstract only: `α ≃ α_loc < α_s`; threshold exists, value unstated |
| **BD / RSOS (1+1) d — Jiang & Xia, CTP 78, 105602 (2026)** | 0–0.4 | "β increases with θ" (values only in figures) | — | — | — | L = 1024; S,K vs θ given; θ_c ≈ 0.2 (distributions), 0.25 (macroscopic) |
| **RD / BD (1+1) d — Wang & Xia, PRE 111, 024124 (2025)** | — | not accessible | not accessible | — | — | paywalled; θ_c = 0.25 reported as emergence threshold |
| **BD-I (1+1) d — Rodríguez-Fernández et al., PRE 110, 024104 (2024)** | 0 → 1/2 | — | increasing with θ | > α at large θ | ≈ 3/2 roughly constant | L = 2^10…2^13; qualitative, figures only |
| **KPZ-II (1+1) d — Rodríguez-Fernández et al. 2024** | < ≈0.23–0.25 | — | 1/2 (plateau) | 1/2 (= α) | — | then α, α_s grow above θ_th ≈ 1/4 |
| **KPZ-I (1+1) d — Rodríguez-Fernández et al. 2024** | all θ > 0 | — | monotonically increasing from θ=0 | > α for large θ | — | no threshold |
| theory: Medina et al. 1989 DRG (1+1) d | < 1/6 | α/z with α=1/2, z=3/2 | 1/2 | — | 3/2 | `α = 1.69θ + 0.22`, `z = (2α+1)/(1+2θ)` for θ>1/6 |
| theory: Squizzato & Canet, PRE 100, 062143 (2019) NPFRG | θ_th | — | — | — | — | **θ_th = 1/6 in d=1; θ_th ≃ 0.346 in d=2** |
| theory: Katzav & Schwartz SCE (2004), Fedorenko FRG (2008) | all θ | — | smooth functions of θ | — | decreasing | **no threshold** |
| experiment/simulation literature (Lam, Sander & Wolf, PRA 46, R6128 (1992)) | θ ∈ [0,1/2) | — | agreed with Medina for some noise generators, not others | — | — | first BD numerics with temporally correlated noise |

---

## 9. PRIOR ART ON HELLINGER DISTANCE AND ISOLATION FOREST IN THIS FIELD

**Hellinger distance — PRIOR ART EXISTS (published, peer-reviewed, and in the identical model class).**
Wang, S. & Xia, H., *Emergence in kinetic roughening with long-range temporal correlations*, **Phys. Rev. E 111,
024124 (2025)**, DOI [10.1103/PhysRevE.111.024124](https://doi.org/10.1103/PhysRevE.111.024124), explicitly states
that "the **Hellinger distance** and the **novelty detection**" are used as "two quantitative measurement methods …
independently to determine whether emergence occurs" for **random deposition and ballistic deposition with
long-range temporal correlations**, and that a threshold exists above which emergence occurs
(threshold = 0.25 per Jiang & Xia 2026). The CUMT institutional press release confirms the Chinese wording
"Hellinger距离和新奇性检测" ([smsp.cumt.edu.cn](https://smsp.cumt.edu.cn/55/6e/c18887a677230/page.htm)).
**A Hellinger-distance-based threshold detection in long-range-temporally-correlated discrete growth is therefore
not novel.**

**Isolation Forest — NO PUBLISHED PRIOR ART FOUND IN THIS FIELD.**
- No paper in this literature that I could access names Isolation Forest. The "novelty detection" of Wang & Xia
  2025 is **not specified** in any accessible source; the full text is paywalled and I could not verify whether it
  is Isolation Forest, One-Class SVM, Local Outlier Factor, or a reconstruction-error method. **This is the single
  most important unresolved question for novelty positioning, and it can only be settled by reading the PRE 111,
  024124 full text (institutional access required — I have none).**
- arXiv full-text search: `"isolation forest" AND "universality"` → 5 hits, all computer science/engineering, none
  physics. `"machine learning" AND "kinetic roughening"` → 0 hits. `"ballistic deposition" AND "machine learning"`
  → 0 hits. `"surface growth" AND "anomaly detection"` → 0 hits.
- The only artefact combining Isolation Forest with surface-growth universality classes is a **public GitHub
  repository** (`adamfbentley/ml-universality-classification`), self-labelled a research prototype and not
  peer-reviewed. If that repository is the present project's own code, it is not citable prior art; if it is
  someone else's, it is still not peer-reviewed.

**Bottom line:** the "distribution-distance metric to locate a universality crossover" idea is **already published
(Hellinger distance, PRE 111, 024124 (2025))**. The distinctive remaining space would have to be (i) *which*
ML anomaly detector is used, (ii) the *feature/embedding* space (raw height fields vs. exponents vs. distributions),
and (iii) a *quantitative comparison of several detectors* — none of which I could confirm as pre-empted, because
the one relevant paper is paywalled.

---

## 10. Sources I could NOT access (explicit list)

| Source | Why | Status |
|---|---|---|
| **PRE 103, 012121 (2021)**, Song & Xia — "Extensive numerical simulations…" | APS paywall; Unpaywall `oa_status: closed`; **no arXiv version exists** | abstract only (PubMed/OpenAlex); **no β(θ), α(θ), θ_th numbers** |
| **PRE 111, 024124 (2025)**, Wang & Xia — "Emergence…" | APS paywall; Unpaywall `oa_status: closed`; no arXiv preprint | abstract only; **novelty-detection algorithm unverified** |
| **J. Stat. Mech. (2016) 113206**, Song & Xia | IOP paywall; abstract elided by publisher in Semantic Scholar; no arXiv version | OpenAlex abstract only |
| Song & Xia JSTAT 2021 — the **numeric θ-dependence of α, β, z** (their Figs. 6, 12) | only plotted, not tabulated | two θ values recovered from text; the rest not readable as numbers |
| Jiang & Xia 2026 — exact `t_1`, `t_sat`, and the saturation time printed as "t = 218" | figure-internal / probable extraction artefact | reported verbatim with caveat |
| sciprofiles.com record for PRE 111, 024124 | HTTP 403 | not used |

**Note on tooling limits (not a source limitation):** `web_fetch` cannot decode `application/pdf`; I worked around
this by downloading PDFs with Python (`urllib`) into the workspace and extracting text with PyMuPDF. IOP's
`iopscience.iop.org/.../pdf` endpoint returns an HTML challenge page to scripted requests; the
`beta.iopscience.iop.org` mirror served the real open-access PDF.

---

## 11. BibTeX

```bibtex
@article{SongXia2021JSTAT,
  author  = {Song, Tianshu and Xia, Hui},
  title   = {Kinetic roughening and nontrivial scaling in the {Kardar--Parisi--Zhang} growth with long-range temporal correlations},
  journal = {Journal of Statistical Mechanics: Theory and Experiment},
  volume  = {2021},
  number  = {7},
  pages   = {073203},
  year    = {2021},
  doi     = {10.1088/1742-5468/ac06c3},
  eprint  = {2102.13356},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.stat-mech},
  note    = {Full text verified (arXiv:2102.13356v2)}
}

@article{SongXia2021PRE,
  author  = {Song, Tianshu and Xia, Hui},
  title   = {Extensive numerical simulations of surface growth with temporally correlated noise},
  journal = {Physical Review E},
  volume  = {103},
  number  = {1},
  pages   = {012121},
  year    = {2021},
  doi     = {10.1103/PhysRevE.103.012121},
  note    = {Paywalled; abstract only verified}
}

@article{SongXia2016,
  author  = {Song, Tianshu and Xia, Hui},
  title   = {Long-range temporal correlations in the {Kardar--Parisi--Zhang} growth: numerical simulations},
  journal = {Journal of Statistical Mechanics: Theory and Experiment},
  volume  = {2016},
  number  = {11},
  pages   = {113206},
  year    = {2016},
  doi     = {10.1088/1742-5468/2016/11/113206},
  note    = {Paywalled; abstract only verified}
}

@article{WangXia2025,
  author  = {Wang, Shuting and Xia, Hui},
  title   = {Emergence in kinetic roughening with long-range temporal correlations},
  journal = {Physical Review E},
  volume  = {111},
  number  = {2},
  pages   = {024124},
  year    = {2025},
  doi     = {10.1103/PhysRevE.111.024124},
  note    = {Paywalled; uses Hellinger distance and novelty detection; algorithm unverified}
}

@article{JiangXia2026,
  author  = {Jiang, Peiyan and Xia, Hui},
  title   = {Height distributions of discrete growth models in the presence of long-range temporal correlations},
  journal = {Communications in Theoretical Physics},
  volume  = {78},
  number  = {10},
  pages   = {105602},
  year    = {2026},
  doi     = {10.1088/1572-9494/ae8f20},
  note    = {Open access CC BY-NC-ND; full text verified}
}

@article{ChangXia2025,
  author  = {Chang, Zhichao and Xia, Hui},
  title   = {Probability distributions for kinetic roughening in the {Kardar--Parisi--Zhang} growth with long-range spatiotemporal correlations},
  journal = {Communications in Theoretical Physics},
  year    = {2025},
  doi     = {10.1088/1572-9494/adfa78},
  eprint  = {2205.08483},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.stat-mech}
}

@article{AlesLopez2019,
  author  = {Al{\'e}s, Alejandro and L{\'o}pez, Juan M.},
  title   = {Faceted patterns and anomalous surface roughening driven by long-range temporally correlated noise},
  journal = {Physical Review E},
  volume  = {99},
  number  = {6},
  pages   = {062139},
  year    = {2019},
  doi     = {10.1103/PhysRevE.99.062139},
  eprint  = {1902.06674},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.stat-mech},
  note    = {Full text verified (arXiv:1902.06674)}
}

@article{RodriguezFernandez2024,
  author  = {Rodr{\'i}guez-Fern{\'a}ndez, Enrique and Al{\'e}s, Alejandro and Mart{\'i}n-{\'A}lvarez, Jorge and L{\'o}pez, Juan M.},
  title   = {Fragility of {Kardar-Parisi-Zhang} universality class in the presence of temporally correlated noise},
  journal = {Physical Review E},
  volume  = {110},
  number  = {2},
  pages   = {024104},
  year    = {2024},
  doi     = {10.1103/PhysRevE.110.024104},
  note    = {Green OA; full text verified via UCrea repository}
}

@article{HuHaoXia2023,
  author  = {Hu, Xiongpeng and Hao, Dapeng and Xia, Hui},
  title   = {Improved finite-difference and pseudospectral schemes for the {Kardar-Parisi-Zhang} equation with long-range temporal correlations},
  journal = {Physica A: Statistical Mechanics and its Applications},
  volume  = {619},
  pages   = {128744},
  year    = {2023},
  eprint  = {2105.13528},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.stat-mech}
}

@article{Medina1989,
  author  = {Medina, Ernesto and Hwa, Terence and Kardar, Mehran and Zhang, Yi-Cheng},
  title   = {Burgers equation with correlated noise: Renormalization-group analysis and applications to directed polymers and interface growth},
  journal = {Physical Review A},
  volume  = {39},
  number  = {6},
  pages   = {3053--3075},
  year    = {1989},
  doi     = {10.1103/PhysRevA.39.3053}
}

@article{SquizzatoCanet2019,
  author  = {Squizzato, Davide and Canet, L{\'e}onie},
  title   = {Kardar-{Parisi}-{Zhang} equation with temporally correlated noise: A nonperturbative renormalization group approach},
  journal = {Physical Review E},
  volume  = {100},
  number  = {6},
  pages   = {062143},
  year    = {2019},
  doi     = {10.1103/PhysRevE.100.062143},
  eprint  = {1907.02256},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.stat-mech},
  note    = {Analytic NPFRG; theta_th = 1/6 (d=1), 0.346 (d=2)}
}

@article{Lam1992,
  author  = {Lam, Chi-Hang and Sander, Leonard M. and Wolf, Dietrich E.},
  title   = {Surface growth with temporally correlated noise},
  journal = {Physical Review A},
  volume  = {46},
  number  = {10},
  pages   = {R6128--R6131},
  year    = {1992},
  doi     = {10.1103/PhysRevA.46.R6128}
}

@article{KatzavSchwartz2004,
  author  = {Katzav, Eytan and Schwartz, Moshe},
  title   = {Kardar-{Parisi}-{Zhang} equation with temporally correlated noise: A self-consistent approach},
  journal = {Physical Review E},
  volume  = {70},
  number  = {1},
  pages   = {011601},
  year    = {2004},
  doi     = {10.1103/PhysRevE.70.011601}
}

@article{Fedorenko2008,
  author  = {Fedorenko, Andrei A.},
  title   = {Elastic systems with correlated disorder: Response to tilt and application to surface growth},
  journal = {Physical Review B},
  volume  = {77},
  number  = {9},
  pages   = {094203},
  year    = {2008},
  doi     = {10.1103/PhysRevB.77.094203}
}
```
