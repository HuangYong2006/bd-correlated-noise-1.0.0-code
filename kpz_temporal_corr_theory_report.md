# KPZ / ballistic deposition with temporally correlated noise — theoretical extraction

## 0. Conventions and access

Noise correlator (Squizzato & Canet, Eq. 4; https://arxiv.org/html/1907.02256v2):

    D(omega,k) = D0 + D_theta * k^(-2 rho) * omega^(-2 theta),        rho = 0 for the pure temporal case

Equivalent real-space form used by the numerical papers (Alés & López, Eq. 3):

    <eta(x,t) eta(x',t')> = 2 D delta_{x,x'} |t - t'|^(2 theta - 1)

d = substrate dimension; alpha = chi (global roughness), beta = chi/z (growth), z (dynamic).
alpha_loc (local) and alpha_s (spectral) are used only by the numerical (Xia-group, Alés-López) papers.

Access: PDFs cannot be fetched by the available tools. arXiv HTML was read via (a) the W3C converter
https://www.w3.org/services/html2txt?url=<encoded> (which can also be chained onto
https://r.jina.ai/ to convert PDFs server-side), and (b) Microlink DOM-node extraction
https://api.microlink.io/?url=<encoded-arxiv-html>&meta=false&data.x.selector=%23S4 (code-block dots
escaped as \.), which bypasses the ~100 kB page cap. Claims I could not verify are labelled explicitly.

---

## 1. Medina, Hwa, Kardar & Zhang (1989), Phys. Rev. A 39, 3053 — one-loop DRG

The paper itself is paywalled and not on arXiv, so everything here is quoted from three independent
papers that reproduce it. (I could not read the original.)

**1.1 Threshold in d = 1: theta_th = 1/6.**
- Squizzato & Canet, Sec. I: "They found that the SR-KPZ fixed point is stable up to a threshold
  value theta_th = 1/6, and thus for theta <= theta_th, the critical exponents are the standard
  SR-KPZ ones z_SR = 3/2 and chi_SR = 1/2." — https://arxiv.org/html/1907.02256v2
- Katzav & Schwartz, Sec. I: "They claim that for phi > 0.167 the correlations become relevant"
  (their phi = theta) — https://arxiv.org/html/cond-mat/0309280v2
- Fedorenko, Sec. V, independently re-derives the same number: "This allows one to compute the exact
  critical value theta_c = 1/6 (for rho = 0) at which there is a crossover from the SR FP to the LR
  FP. For arbitrary rho, the crossover to the LR FP happens at 6 theta + 4 rho > 1. (77)"
  — https://arxiv.org/html/0712.0801v2

**1.2 One-loop long-range exponents — TWO equivalent printed forms.**
- Squizzato & Canet, Eq. (5): chi_LR = (1 + 4 theta)/(3 + 2 theta), z_LR = 2 - chi_LR
  — https://arxiv.org/html/1907.02256v2
- Fedorenko, Eqs. (79)-(80), quoting the same DRG result for general rho:
  z*(theta,rho) = 2 - (1 + 4 theta + 2 rho)/(3 + 2 theta), chi*(theta,rho) = (1 + 4 theta + 2 rho)/(3 + 2 theta);
  at rho = 0 this is z* = 5/(3 + 2 theta), chi* = (1 + 4 theta)/(3 + 2 theta)
  — https://arxiv.org/html/0712.0801v2
  (My algebra: z = 2 - chi and z = 5/(3+2 theta) are identical.)
- Validity caveat, verbatim: the one-loop expression was "obtained by neglecting the corrections on
  the non-linearity induced by the violation of Galilean invariance due to the temporal correlations.
  This expression is thus only valid for small theta close to the threshold."
  — https://arxiv.org/html/1907.02256v2

**1.3 Numerical solution of the truncated one-loop flow equations (d = 1), fitted:**
    chi_LR = 1.69 theta + 0.22,     z_LR = (2 chi_LR + 1)/(1 + 2 theta)
- Squizzato & Canet Eq. (7): https://arxiv.org/html/1907.02256v2
- Katzav & Schwartz Eqs. (6)-(7) print the same two formulas as alpha_DRG(phi) = 1.69 phi + 0.22 and
  z_DRG(phi) = (2 alpha_DRG(phi) + 1)/(1 + 2 phi), also attributing them to Medina et al.
  — https://arxiv.org/html/cond-mat/0309280v2
- Fedorenko Eq. (83) calls this the "modified DRG" result and shows it violates his stability bound
  (see §4): "the exponent (83) computed using the modified DRG violates the criterion of the LR FP
  stability, and thus is ruled out" — https://arxiv.org/html/0712.0801v2
- Katzav & Schwartz add: the DRG approach "found no solution above phi = 0.46 (claiming that no stable
  surface can grow under the condition 0.46 < phi < 0.5)" — https://arxiv.org/html/cond-mat/0309280v2

**1.4 The exact long-range scaling relation.**
    z_LR (1 + 2 theta) - 2 chi_LR - d = 0                       (Squizzato & Canet, Eq. 6)
quoted as "the exact relation" that replaces the Galilean relation chi_SR + z_SR = 2 at the LR fixed
point — https://arxiv.org/html/1907.02256v2. Fedorenko's version, for general rho, is
    z_LR (1 + 2 theta) - 2 chi_LR + 2 rho = 1                    (his Eq. 78)
and it follows from his Eq. (74), d_l U_theta = U_theta [ z (1 + 2 theta) - 2 chi - 1 + 2 rho ], which
he states is EXACT TO ALL ORDERS because U_theta acquires no fluctuation corrections
— https://arxiv.org/html/0712.0801v2. See §6 for the power-counting derivation.

**1.5 One sentence in Squizzato & Canet that appears to be an error.** Their text says the exact
relation (6) "is violated by the estimate (5)". Substituting (5) into (6) gives zero identically
(verified algebraically):
    (2 - chi)(1 + 2 theta) - 2 chi - 1 = 1 + 4 theta - chi(3 + 2 theta) = 0.
The one-loop estimate as printed DOES satisfy the exact LR relation for every theta (numerically
checked at theta = 1/6, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5; residual exactly 0). Fedorenko's
independent printing of the same one-loop result in the form z* = 5/(3+2 theta), chi* = (1+4 theta)/(3+2 theta)
confirms this. Do not use that sentence to argue that the one-loop exponents are internally inconsistent.

**1.6 A second special value of theta: 1/4 (a DRG singularity, not the threshold).**
- Squizzato & Canet, Sec. IV.2.2: "In the work of [18] (Medina et al.), a non-physical divergence in
  omega appears for theta >= 1/4. Within the present work, a similar singularity at theta = 1/4 arises
  in the limit of zero frequency ... in the flow equation of f^D_kappa, in the term proportional to
  (omega + varpi)^(-4 theta)." — https://arxiv.org/html/1907.02256v2
- Alés & López, Discussion: "the value theta = 1/4 was already shown to play a special role in the
  perturbative RG approximation of Medina et al., as the point at which the renormalized noise
  amplitude D*(omega) has a singular correction at leading omega-order. Further singularities appear at
  larger values of theta: making the RG treatment ill-constructed."
  — https://ar5iv.labs.arxiv.org/html/1902.06674
- Fedorenko also notes (quoting Medina et al.) that the LR fixed point in the truncated flow exists
  only for theta < 1/4 — https://arxiv.org/html/0712.0801v2

---

## 2. Squizzato & Canet (2019), Phys. Rev. E 100, 062143; arXiv:1907.02256 — NPRG

All quotes read directly from the paper's arXiv HTML: https://arxiv.org/html/1907.02256v2
(ar5iv mirror: https://ar5iv.labs.arxiv.org/html/1907.02256)

**2.1 Thresholds.**
- d = 1: theta_th = 1/6. Sec. IV.2.1: "For theta < theta_th = 1/6, g_kappa flows to a finite fixed
  point value g_* while w^theta_kappa flows to zero. At the same time, eta^lambda_kappa also flows to
  zero, hence Galilean invariance is dynamically restored, and the critical exponents take the SR-KPZ
  values chi_SR = 1/2 and z_SR = 3/2." ... "The critical value theta_th = 1/6 delimiting the two
  regimes is clearly identified on Fig. 3 by the algebraic decay of w^theta_kappa and eta^lambda_kappa
  with the RG time s." Fig. 3 caption: "The critical value is theta_th = 0.166 confirming the
  theoretical prediction theta_th = 1/6."
  For theta > theta_th: "both g_kappa and w^theta_kappa flow to a non-zero fixed-point value, and the
  violation of Galilean invariance eta^lambda_* increases with theta ... Hence in this regime, the
  long-distance properties are controlled by a line of LR fixed points, with critical exponents
  depending on theta."
- d = 2: theta_th = 0.346. Sec. IV.2.2: "Substituting the NLO value of eta^nu_* = 2 - z = chi in (49),
  one obtains a theoretical estimate of the critical value theta_th = 0.346." (Introduction says
  "theta_th ~= 0.35 in d = 2".)
- Their finding is stated to confirm the earlier RG results: "These findings are in agreement with the
  results presented in [18, 37], and show that a pure SR-KPZ regime is not destroyed for an
  infinitesimal theta, contrary to the scenario advocated by SCE or Flory approaches."
  ([18] = Medina et al., [37] = Fedorenko.)

**2.2 NUMERICAL VALUES — what actually exists in the paper.**
There is NO table of exponents in Section IV, and the paper NEVER reports the growth exponent beta
(verified: zero occurrences of "Table", "beta", "alpha" in Sections I-IV). chi(theta) and z(theta) are
given as Figures 4 (d = 1) and 6 (d = 2). The explicit numbers are:

- d = 1, chi in the LR phase (NLO), Fig. 4 caption verbatim: "At NLO, we find an almost linear
  behavior of chi(theta) ~= 1.43 theta + 0.26. It turns out to be reasonably close to numerical
  estimation from the approximate RG equations (7) for chi, but not for z. All the estimates from NPRG
  lie within the bounds given in (9) (shaded region)." (Their eq. (7) = Medina's fit 1.69 theta + 0.22.)
- d = 1, z in the LR phase, Sec. IV.2.1 verbatim: "[more quantitatively, z_LR = 1.42(5) for theta = 0.47
  and z_LR = 1.49(4) for theta = 0.24]". These are their own Fig. 4(b) NPRG values (Fig. 4 contains
  only their NLO / NLO^omega / LPA curves plus the DRG one-loop curve; no simulation points).
- d = 1, Sec. IV.2.1: "Let us emphasize that all our estimates for z_LR(theta) are decreasing with theta
  and lie within the bounds (9) derived in [37]" (bounds (9) = Fedorenko's, see §4).
- Summary of what is printed for d = 1: (theta <= 1/6: chi = 1/2, z = 3/2); (theta = 0.24: z = 1.49(4));
  (theta = 0.47: z = 1.42(5)); plus the chi fit above. Between theta_th and 0.47 the curves are graphical only.
- d = 2: "the estimation of chi for the pure SR-KPZ within the NLO approximation is chi_SR ~= 0.375, in
  good agreement with the value chi ~= 0.3869(4) from recent numerical simulations [72]" (0.3869(4) is
  the SIMULATION value). They use the SR relation eta^nu_* = 2 - z = chi, so z ~= 1.625 at their SR
  point (my arithmetic). LR exponents in d = 2 exist ONLY as curves in Fig. 6 (two NLO variants, "Both
  results are in agreement"); no numbers, no fit, no table.
- Consistency with the exact LR relation (my arithmetic, labelled as such): inserting their chi fit
  chi = 1.43 theta + 0.26 into z = (2 chi + 1)/(1 + 2 theta) (Eq. 6 with d = 1) gives z = 1.4908 at
  theta = 0.24 and z = 1.476 at theta = 0.47 — consistent with their printed 1.49(4) and 1.42(5) within
  errors. The same fit gives chi = 0.4983 ~= 1/2 at theta = 1/6, i.e. it passes through the SR point.
- Robustness checks they report: LPA gives the same theta_th with exponents closer to the one-loop
  results; the NLO^omega scheme has a residual Galilean violation at theta = 0,
  eta^lambda_* ~= 0.0065, which slightly shifts theta_th (other exponents stay close to NLO).
- Relation they use in d = 1 for its own scheme: "z = 2 - chi - eta^lambda_*" (the Galilean relation
  chi + z = 2 "only holds at the SR fixed point").

**2.3 Anomalous scaling (Sec. IV.3, "Anomalous scaling at the LR fixed point").**
- They use the term "intrinsic anomalous scaling" and never use alpha_s or alpha_loc. Verbatim:
  "This anomalous scaling corresponds to what is termed intrinsic anomalous scaling in [74], for which
  chi < 1 but a new independent exponent exists, by opposition to superroughnening, for which chi > 1."
- Two exponents: the usual chi^(kappa) and chi^(p) obtained "in which the scaling dimensions eta^X_*
  are replaced by the actual scaling exponents in momentum space", with
  chi^(p) = (2 - d + eta^D_{p>>1} - eta^nu_{p>>1})/2, and
  zeta^D(y) = zeta^D_0 y^(-Delta eta^D) for y -> 0, zeta^D_inf y^(-eta^D_*/z) for y -> infinity,
  Delta eta^D = eta^D_* - eta^D_{p>>1}.
- "Hence the quantity Delta chi = chi^(kappa) - chi^(p) is a direct measurement of the presence of
  anomalous scaling in the system."
- "Although in our case, the anomalous scaling appears in the whole LR phase theta > theta_th, it
  becomes significant only at larger values of theta >~ 0.3, which is consistent with the result of
  [41] which reports large anomalous scaling for theta >~ 0.3 - 0.4. However, for such large values,
  the numerical values for chi in [41] are greater than one, which would correspond to superroughening,
  whereas our results remains in the regime chi < 1, which corresponds to intrinsic anomalous scaling."
- "The value of eta^D_{p>>1} cannot be predicted by any scaling argument."
- In the LR phase "z + chi != 2"; and the new exponent is linked to "emergent spatial structures in the
  system, referred to as faceting, which appears for theta >~ 0.23" (0.23 = Alés & López's threshold).

**2.4 Their explicit statement of disagreement with the simulations (Sec. IV.2.1), verbatim:**
  "Let us finally note that these theoretical results are not in agreement with the most recent numerical
  simulations [41]. In the simulations, the SR-KPZ phase is destroyed for any non-zero theta (no threshold
  value) and the value of the dynamical critical exponent is found to be constant z ~= 3/2 independently
  of theta. We have no clear understanding of these discrepancies. The main differences lie in the fact
  that the simulations deal with finite-size systems in discretized space and time whereas we work in the
  infinite-size limit in continuous space-time. In general, finite-size effects tend to smear out
  transitions. A putative explanation could be that very large system sizes are required to clearly
  resolve the threshold." They add that z depends only on eta^nu_*, "which is less sensitive than
  eta^D_* to the presence of the LR noise", so the numerical z ~= 3/2 may still be compatible.

**2.5 Canonical dimensions they state (Eq. 12), verbatim:** "[g] = 2 - d, [w_theta] = 4 theta", with
g = lambda^2 D0/nu^3 and w_theta = D_theta/(D0 nu^(2 theta)).

---

## 3. The "no threshold" alternatives

**3.1 Ma & Ma (1993), Flory-type scaling analysis.**
Citation: Hanfei and Ben-Kun Ma, "Scaling analysis of Langevin-type equations", Phys. Rev. E 47, 3738
(1993), DOI 10.1103/PhysRevE.47.3738 — cited in the literature as "Ma & Ma (1993)" and, by Jiang & Xia,
as "Han & Ma"; the first author appears as the single name "Hanfei" in the ADS/Semantic Scholar record
(https://api.semanticscholar.org/graph/v1/paper/DOI:10.1103/PhysRevE.47.3738?fields=title,authors,year,externalIds).
Their exponents as reproduced by two independent papers:

    chi_LR = (2 + 4 theta)/(2 theta + d + 3),      z_LR = (2d + 4)/(d + 3 + 2 theta)

- Squizzato & Canet, Eq. (8), with the comment "such that the SR-KPZ exponents are only recovered at
  theta = 0" — https://arxiv.org/html/1907.02256v2
- Katzav & Schwartz, Eqs. (8)-(9), as alpha_SA(phi) = (2 + 4 phi)/(2 phi + d + 3) and
  z_SA(phi) = (2d + 4)/(2 phi + d + 3) — https://arxiv.org/html/cond-mat/0309280v2
In d = 1 (my algebra): chi = (1 + 2 theta)/(2 + theta), z = 3/(2 + theta), beta = (1 + 2 theta)/3.
Two checks (my algebra, numeric): this branch (i) violates the exact LR relation
z(1 + 2 theta) = 2 chi + d for every theta < 1 (residual -0.385 at theta = 1/6, -0.20 at theta = 1/2),
and (ii) lies BELOW Fedorenko's lower bound z >= 5/(3 + 2 theta) for every 0 < theta < 1 (crossing only
at theta = 1). Squizzato & Canet state that the bound "rules out both the second SCE solution and the
scaling solution" (the "scaling solution" = these Flory exponents).

**3.2 Katzav & Schwartz (2004), self-consistent expansion (SCE).**
Citation: Phys. Rev. E 70, 011601 (2004), arXiv:cond-mat/0309280. Source for all quotes:
https://arxiv.org/html/cond-mat/0309280v2
Their notation: alpha = (Gamma - d)/2, z(Gamma) = (d + 4 - Gamma)/2, beta = alpha/z; noise
<eta(r,t) eta(r',t')> = 2 D0 delta(r - r') |t - t'|^(2 phi - 1) with phi < 1/2.

They split the (Gamma, z) quadrant into sectors and find, in the relevant sector ("delta"), TWO
strong-coupling solutions:
1. a transcendental solution from their Eq. (60), denoted Gamma_phi(d), z_phi(d), which "reduces to the
   standard KPZ results when phi = 0". It exists only if Gamma > z(1 + 2 phi), i.e.
   phi < [3 Gamma_phi(d) - d - 4] / [2(d + 4 - Gamma_phi(d))];
2. an EXPLICIT solution, valid for phi larger than that bound:

       z = (d + 4)/(3 + 2 phi),        Gamma = (d + 4)(1 + 2 phi)/(3 + 2 phi),

   i.e. in d = 1 (using alpha = (Gamma - 1)/2):

       alpha = (1 + 4 phi)/(3 + 2 phi),   z = 5/(3 + 2 phi),   beta = alpha/z = (1 + 4 phi)/5.

   (Note: this explicit SCE branch is numerically identical to Medina et al.'s one-loop formula of §1.2.)
They conclude (verbatim): "the actual strong coupling exponents of the KPZ problem with temporally
correlated noise are just Gamma = max{Gamma_phi(d), (d + 4)(1 + 2 phi)/(3 + 2 phi)} and its
corresponding z. Thus, the transition between the two solutions as a function of phi (if such a
transition exists) is continuous. However, it should be emphasized that for a specific phi one of these
solutions dominates so there is no phase transition between them."
Numbers they print: Gamma_0(1) = 2 analytically at phi = 0 in d = 1 (hence alpha_0 = 1/2, z_0 = 3/2);
Gamma_0(2) = 2.59 numerically in d = 2; the crossover value phi_c = [3 Gamma_0(d) - d - 4]/[2(d + 4 - Gamma_0(d))]
equals 1/6 for d = 1. There is NO table of exponents vs phi; the phi-dependence appears only in their
Fig. 2 (alpha_phi, beta_phi, z_phi at d = 1).
Key qualitative claims, verbatim:
- "using the self-consistent approach we found no 'threshold behavior'. That is, we found a continuous
  variation of the scaling exponents alpha and z as a function of phi over the whole range of possible
  phi's, rather than no variation of these exponents up to a critical value of phi_c and a quasi-linear
  behavior from that point on."
- "Third, we found that z is a decreasing function of phi, while the DRG approach predicts an increasing
  value of z."
- On the coincidence of thresholds: "Interestingly, the threshold phi_c (the crossover point) that was
  predicted using DRG in d = 1 (namely 0.167) is the same as the lower bound we found above (phi_c = 1/6
  for d = 1). Our exact statement was that for phi > 1/6 the second strong-coupling KPZ solution becomes
  possible in principle (but not in practice). Therefore, the DRG result might reflect this exact
  statement." and "Generally speaking, the two methods disagree on the values of the critical exponents
  significantly over most of the phi-range."
- In d = 1 they also write that the explicit branch is "practically irrelevant, since Gamma_phi dominates
  the whole phi range", but in a later paragraph state "the actual solution that determines the exponents
  is the second solution (in one dimension)". These two statements appear to contradict each other; I
  report both. Squizzato & Canet's reading is that the SCE's DOMINANT solution is the one "which leads to
  a smooth dependence on theta with no threshold, and with a decreasing z_LR(theta), whereas the solution
  (7) [Medina's fit] is increasing" — https://arxiv.org/html/1907.02256v2. Either way, the dominant SCE
  branch is the transcendental one, for which no closed form exists (numerical solution only).
- Fedorenko reports that his stability bound rules out "the second SCA solution considered to be
  dominant" — https://arxiv.org/html/0712.0801v2.

---

## 4. Fedorenko (2008), Phys. Rev. B 77, 094203; arXiv:0712.0801 — two-loop results and the bound

Source: https://arxiv.org/html/0712.0801v2 (full text obtained; journal ref and arXiv id verified at
https://arxiv.org/abs/0712.0801 and
https://api.semanticscholar.org/graph/v1/paper/DOI:10.1103/PhysRevB.77.094203?fields=title,externalIds,year).

**4.1 Framework.** Double expansion in epsilon = 4 - d and delta = 4 - a (a = LR disorder correlation
exponent). The KPZ interface maps to a directed polymer, "a one-dimensional (d = 1, epsilon = 3) elastic
object", with the temporal-noise exponent related by "relating delta = 3 + 2 theta".

**4.2 Two-loop ingredients (as printed).**
- psi = -(1/4)(delta - epsilon) R_2^(4)(0)  [Eq. 30]; since R_2^(4)(0) > 0, psi < 0 ("psi is nonpositive
  for columnar disorder").
- zeta_LRRB = (delta + 2 psi)/5  [Eq. 34, exact identity]; zeta_LRRF = (delta + 2 psi)/3 [Eq. 37];
  one-loop: zeta_LRRB = delta/5.
- Improved two-loop psi: psi = epsilon^2 f1(delta/epsilon) + epsilon^3 f2(delta/epsilon) + ... [Eq. 35]
  with f2(x) = [(0.159192 x - 0.200087)/(x - 1.04149)] f1(x) [Eq. 36].
- Mapping to the dynamic exponent: z(d') = 1/zeta(d = 1, N = d') [Eq. 86], giving

      z_LR = (5 - 2 rho)/(3 + 2 theta + 2 psi)                                    [Eq. 87]

  CAVEAT: the plain-text renderings of Eq. (87) display the denominator as "3 + 2 theta + psi"
  (the coefficient of psi is lost in text extraction). The 2 psi form is the one consistent with
  (i) reproducing Eq. (77) exactly when psi = 0, z = 3/2, (ii) eliminating psi between Eqs. (87) and
  (90) to reproduce Eq. (78), and (iii) zeta_LRRB = (delta + 2 psi)/5 with Eq. (86). I flag this as
  a text-extraction ambiguity; check the printed PRE page if you need the exact coefficient.
- Other exact relations printed: Eq. (74) d_l U_theta = U_theta [z(1 + 2 theta) - 2 chi - 1 + 2 rho]
  (exact to all orders because U_theta is not renormalized); Eq. (78)
  z_LR(1 + 2 theta) - 2 chi_LR + 2 rho = 1; Eq. (90) chi + z = 2 - psi z.
- One-loop/DRG comparison values printed as Eqs. (79)-(80): z* = 2 - (1 + 4 theta + 2 rho)/(3 + 2 theta),
  chi* = (1 + 4 theta + 2 rho)/(3 + 2 theta). "To one-loop order in FRG, i.e., for psi = 0, exponent (87)
  coincides with the estimation given by DRG (79) for small theta."

**4.3 The bound and the threshold (Sec. V).**
- Bound, verbatim: "Taking into account that psi is nonpositive for columnar disorder, we obtain the
  lower and upper bounds on z(theta) for rho = 0 and theta in [1/6, 1/2] as
      5/(3 + 2 theta) <= z(theta) <= 3/2.   (88)"
- Derivation, verbatim: "The stability criterion assures that the LR FP in the FRG picture is stable if
  zeta_LR > zeta_SR. This implies that the noise temporal correlations in surface growth problem are
  relevant only if the corresponding dynamic exponents fulfill the condition z_LR < z_SR. Note that this
  criterion is purely based on the mapping between the DP and KPZ problems."
- Threshold, verbatim: "This allows one to compute the exact critical value theta_c = 1/6 (for rho = 0)
  at which there is a crossover from the SR FP to the LR FP. For arbitrary rho, the crossover to the LR
  FP happens at 6 theta + 4 rho > 1. (77)" He attributes the exactness of Eq. (74) to Medina et al.
  (his ref. 52) and re-derives the crossover from z_LR < z_SR = 3/2.
- What it rules out, in his own words: (a) Medina's "modified DRG" fit of §1.3 — "Since
  z_SR(d'=1, rho=0) = 3/2 the exponent (83) computed using the modified DRG violates the criterion of
  the LR FP stability, and thus is ruled out"; (b) "the second SCA solution considered to be dominant
  does not satisfy bounds (88), and thus is ruled out". Squizzato & Canet add that the bound also
  excludes the Flory/"scaling" solution (Ma & Ma), which Fedorenko's own text does not name.
- NUMBERS AT epsilon = 3 (d = 1): none. He writes "the large value epsilon = 3 of the expansion
  parameter describing the DP problem makes the estimation of psi highly unreliable", and the Fig. 8
  caption repeats that epsilon = 3 "does not allow for accurate computation of psi and z at large theta".
  z(theta) at eps = 3 exists only as curves in Fig. 8 (not text-extractable). The only eps = 3 number
  printed is for the short-range RB class: zeta_SRRB(eps = 3) = 0.6866 vs the exact 2/3 (Eq. 31), plus an
  improved Pade form (Eq. 32). There is no closed-form two-loop epsilon-expansion of z_LR(theta): psi is
  a numerical fixed-point quantity.

---

## 5. NEGATIVE theta (noise suppressed at low frequency)

No paper I could access treats theta < 0 as its main subject; all parametric studies restrict to
theta >= 0 (Squizzato & Canet's figures use theta = 0 ... 0.47; Alés & López: "The exponent theta in
[0, 1/2) characterizes the temporal correlation range of the noise that becomes more long-term
correlated as theta is increased from zero", https://ar5iv.labs.arxiv.org/html/1902.06674; Katzav &
Schwartz assume phi < 1/2, https://arxiv.org/html/cond-mat/0309280v2). Fedorenko 2008 contains no
statement about theta < 0 (full-text search, https://arxiv.org/html/0712.0801v2).

What IS verifiable:
1. Katzav & Schwartz, Sec. V (their sector "alpha"), verbatim: "...the condition 2 - Gamma < -2 z phi can
   be met only for phi <= 0. Namely, for the case of noise anticorrelated in time [27]. Such a solution
   holds only for d > 4(1 - 2 phi)." and "Combining the results for sectors alpha and beta, we see that
   for phi <= 0 the noise term is irrelevant, and the critical exponents that describe the EW problem
   with uncorrelated noise (namely Gamma = z = 2) are restored. This option is possible if d > 2. In
   addition, for phi > 0 we get the new solution z = 2 and Gamma = 2 + 4 phi ... only for d > 2 + 4 phi.
   Therefore, the lower critical dimension in this problem is d_c = 2 + 4 phi (provided phi > 0,
   otherwise d_c = 2 as mentioned above)."
   URL: https://arxiv.org/html/cond-mat/0309280v2
   Reading: for phi <= 0 the temporally correlated part of the noise is IRRELEVANT and one recovers the
   uncorrelated-noise exponents — in d > 2 that is the EW (linear) fixed point; in d <= 2 the EW fixed
   point is unstable, so the statement translates into SR-KPZ.
2. Squizzato & Canet Eq. (12): "the LR correlation introduces another dimensionless coupling
   w_theta = D_theta/(D0 nu^(2 theta)) ... [g] = 2 - d, [w_theta] = 4 theta"
   (https://arxiv.org/html/1907.02256v2). With their convention (positive canonical dimension =
   relevant, as for g = lambda^2 D0/nu^3 whose dimension 2 - d makes it relevant for d < 2), w_theta is
   relevant for theta > 0 and IRRELEVANT for theta < 0.
3. Deduction for the d = 1 case you care about (mine, not a quotation): for theta < 0
   D(omega) = D0 + D_theta omega^(2|theta|) is analytic at omega = 0, and the real-space correlator
   |t - t'|^(2 theta - 1) decays faster than 1/|t - t'|, i.e. is short-ranged in time; therefore SR-KPZ
   must control the long-time behaviour: chi = 1/2, z = 3/2, beta = 1/3 in d = 1, with no crossover.
   Seen from Fedorenko's Eq. (87)/(77) with theta < 0, the LR branch z_LR = 5/(3 + 2 theta) exceeds
   z_SR = 3/2, i.e. it fails the relevance criterion z_LR < z_SR, so again no LR phase. The only
   practical caveat is UV: D_theta omega^(2|theta|) grows at large omega, which lattice simulations
   regularize automatically.

---

## 6. Power counting: where 1/6 comes from (and where it does not)

**6.1 Naive Gaussian power counting gives NO threshold.** With x -> b x, t -> b^z0 t, h -> b^chi0 h and
the Gaussian (EW) values chi0 = (2 - d)/2, z0 = 2, the spectral noise density has dimension
[D] = 2 chi0 - z0 + d = 0, and the LR part gives

    [D_theta] = 2 chi0 - z0 (1 + 2 theta) + d = (2 - d) - 2(1 + 2 theta) + d = -4 theta,

i.e. |[w_theta]| = 4 theta in the convention of Squizzato & Canet, Eq. (12) ("[g] = 2 - d,
[w_theta] = 4 theta", https://arxiv.org/html/1907.02256v2). In that convention positive dimension =
relevant, so the temporal correlation is canonically RELEVANT for any theta > 0 — naive power counting
predicts no threshold (the same conclusion as the Flory and SCE scenarios of §3).

**6.2 The exact LR scaling relation (three lines).** Requiring the LR noise term to be marginal at the LR
fixed point, [D_theta] = 0 with the actual exponents chi, z, gives

    2 chi - z (1 + 2 theta) + d = 0      <=>      z (1 + 2 theta) = 2 chi + d,

which is exactly Squizzato & Canet's Eq. (6) and, for general rho, Fedorenko's Eq. (78) (obtained from
his Eq. (74), which is exact to all orders). URLs: https://arxiv.org/html/1907.02256v2 ,
https://arxiv.org/html/0712.0801v2

**6.3 The threshold theta_th = 1/6, now with a citable derivation.** Fedorenko's stability criterion is
"the noise temporal correlations ... are relevant only if ... z_LR < z_SR", which for rho = 0 gives the
crossover line 6 theta + 4 rho = 1, i.e. theta_c = 1/6
(https://arxiv.org/html/0712.0801v2, Eq. (77)). Equivalently, in the language of §6.2, evaluating the
marginality condition at the SR-KPZ fixed point (chi_SR = 1/2, z_SR = 3/2 in d = 1) gives

    z_SR (1 + 2 theta) = 2 chi_SR + d   =>   (3/2)(1 + 2 theta) = 2   =>   theta = 1/6,

which is also exactly where Medina's LR branch chi_LR = (1 + 4 theta)/(3 + 2 theta) meets chi = 1/2
(and z = 5/(3 + 2 theta) meets 3/2, and Fedorenko's lower bound meets his upper bound). All three
statements are elementary algebra on formulas quoted in §1 and §4; they are consistent with each other,
and the resulting 1/6 is the threshold reported by Medina et al., by Fedorenko (independently
re-derived) and by Squizzato & Canet's NPRG flow.

---

## 7. Continuum problem vs the DISCRETE ballistic deposition model

**7.1 The discrete model actually simulated.** Alés & López, "Models"
(https://ar5iv.labs.arxiv.org/html/1902.06674): ballistic deposition with

    h(x, t+1) = Max[ h(x,t) + eta(x,t), h(x-1,t), h(x+1,t) ],

integer heights, eta in {0,1}, periodic boundaries, parallel update (all even sites at even times).
The deposition sequence at each site is obtained by thresholding a Mandelbrot fast-fractional-Gaussian-
noise sequence: "In the case of the simulations of particle deposition by ballistic deposition we map
eta(x,t) = 1 if Z_x(t) > 0 and eta(x,t) = 0 otherwise. This digitalization of the noise enhances the
statistics of the simulations by avoiding the formation of overhangs on the surface." The noise is
spatially uncorrelated with a power-law temporal correlator <eta(x,t) eta(x',t')> = 2D delta_{x,x'}
|t - t'|^(2 theta - 1) and power spectrum <|eta_hat(x,omega)|^2> ~ omega^(-2 theta) for omega -> 0,
theta in [0, 1/2). Jiang & Xia use the same generator family (FFGN of Mandelbrot 1971 plus the binary
mapping of Lam et al. 1992; parameters B = 2, a = 6, u_n = a B^(-n), r_n = exp(-u_n), N = 40 components;
"the FFGN method becomes numerically unstable for large temporal correlation exponents theta > 0.4"),
with L = 1024 and both BD and RSOS rules (BD: particles stick on first contact; RSOS: deposition
accepted only if all height differences |Delta h| <= 1).
URLs: https://ar5iv.labs.arxiv.org/html/1902.06674 , https://beta.iopscience.iop.org/article/10.1088/1572-9494/ae8f20/meta

**7.2 Discrete-model results (Alés & López 2019, PRE 99, 062139; arXiv:1902.06674).**
- "Figure 4 clearly shows that the roughness exponents split up at theta_th = 0.25 +- 0.03 for the
  ballistic deposition model (a similar threshold value theta_th ~= 0.23 was found for the KPZ equation)."
- theta = 0.15: standard scaling, alpha = alpha_s = 0.56 +- 0.02, z = 1.50 (data-collapse value).
- theta = 0.47 (faceted regime): alpha = 1.02 +- 0.01, alpha_s = 1.23 +- 0.02, z = 1.50.
- Columnar limit theta = 1/2 (values from the cited Szendro et al. 2007): alpha = 1.07 +- 0.05,
  alpha_s = 1.5 +- 0.05.
- Comparison with theory: "while the global roughness exponent is nicely predicted by the dynamic RG
  calculations of Medina et al. above the threshold, it fails to describe the correlation effects below
  theta_th. Obviously, no present theory is capable of predicting the existence of facets and the
  associated spectral roughness exponent as the noise correlation range is increased."
- Their Fig. 4 compares with three theory curves: Medina's DRG, the Flory/SA of Hanfei & Ma (1993), and
  Katzav & Schwartz's SCE.
- Mechanism: faceting from exponential localization of phi = exp h around random centers (the theta = 1/2
  columnar limit, per Szendro et al.), extended to theta < 1/2; "the infinitely many singularities [of
  the RG] are a reflection of the failure of the RG treatment to describe the appearance of a new
  exponent alpha_s != alpha."
- They state that the anomalous spectral exponent alpha_s leaves no trace in the usual height-height
  correlation function or the width W(L,t) — it is visible only in the structure factor.

**7.3 Earlier discrete simulations.** Katzav & Schwartz (https://arxiv.org/html/cond-mat/0309280v2):
"The DRG predictions have been checked numerically by Lam et al. [17] using the Ballistic Deposition
model. They found sensible agreement between the DRG prediction and the numerical values they obtained.
However, substantial deviations were found, centered around the expected threshold point phi_0 = 0.167.
Thus, the authors believe that these discrepancies are due to a crossover effect in the simulation and
not due to any approximation in the DRG calculation." [17] = C.-H. Lam, L. M. Sander, D. E. Wolf,
"Surface growth with temporally correlated noise", Phys. Rev. A 46, R6128 (1992) (title/DOI verified at
https://api.semanticscholar.org/graph/v1/paper/DOI:10.1103/PhysRevA.46.R6128?fields=title,authors,year).
They also note that the SCE and the simulations show no threshold either.

**7.4 Recent discrete-model work.** P.-Y. Jiang and H. Xia, "Height distributions of discrete growth
models in the presence of long-range temporal correlations", Commun. Theor. Phys. 78, 105602 (2026),
DOI 10.1088/1572-9494/ae8f20 (full text at
https://beta.iopscience.iop.org/article/10.1088/1572-9494/ae8f20/meta). BD and RSOS, theta in [0, 0.4]:
- "We find a critical threshold theta_c = 0.2 in discrete growth models, consistent with observations for
  the continuous KPZ equation [58,59]. However, this value differs from the threshold theta = 0.25
  derived from macroscopic dynamic phenomena [46,60], suggesting that height fluctuations are more
  sensitive to long-range temporal correlations."
- Height distributions: TW-GOE-like for 0 <= theta <= 0.2 (S = 0.2934, K = 3.1652), with a pronounced
  increase of S and K for theta > 0.2 (RSOS: increase over 0.2 <= theta <= 0.4); saturation
  height-increment distributions converge to Baik-Rains (S = 0.3594, K = 3.2891).
- They report NO exponent values: alpha, alpha_loc, alpha_s and z do not appear anywhere in the paper,
  and beta appears only as "the growth exponent beta is found to increase with theta" (figure inset).
  There are no tables in the paper.
- They also state "this coupling breaks Galilean invariance within the KPZ framework, manifesting as a
  violation of the scaling identity alpha + z = 2".
- Theory citations in that paper attach no numbers to Medina (they quote "theta_c = 0.167"), to
  Hanfei & Ma, or to Katzav & Schwartz; Squizzato & Canet appear only inside a range citation [43-47];
  Fedorenko is not cited.

**7.5 Continuum vs discrete may not share the same theta.** X. Hu, D. Hao and H. Xia, Physica A 619,
128744 (2023) = arXiv:2105.13528, abstract verbatim: "Furthermore, our investigations suggest that when
the effects of long-range temporal correlation are present, the continuum and discrete growth systems do
not belong to the same universality class with the same temporal correlation exponent."
URL: https://arxiv.org/html/2105.13528v2
The same paper is the one that actually uses the alpha / alpha_loc / alpha_s framework; at theta = 0.30
in 1+1 dimensions they report (their own simulations): "we obtain the growth exponents beta = 0.583 +- 0.008
for the PS scheme, and beta = 0.580 +- 0.010 for LS case ... we obtain the global roughness exponents
alpha = 0.73 +- 0.03 for PS and alpha = 0.73 +- 0.03 for LS scheme", with alpha_loc = 0.788(5) / 0.783(5)
and alpha_s = 0.828(6) / 0.827(6). Their abstract states the general pattern: "three characteristic
roughness exponents (global roughness exponent alpha, local roughness exponent alpha_loc, and spectral
roughness exponent alpha_s) are approximately equal within the small temporally correlated regime, and
satisfy alpha_loc ~= alpha < alpha_s for the large temporally correlated regime."
URL: https://ar5iv.labs.arxiv.org/html/2105.13528

**7.6 Discrete-model-related numbers in d = 2 (for completeness).** T.-S. Song and H. Xia, J. Stat. Mech.
2016, 113206 (closed access, no arXiv version). Their values are quoted in the Xia-group follow-up
(arXiv:2306.06952): 2+1-dimensional KPZ with (a) Gaussian noise: alpha = 0.385, z = 1.571; (b) long-range
temporal correlations theta = 0.30: alpha = 0.669, z = 1.474.

**7.7 Summary of the discrete-vs-continuum tension.**
- The DRG/NPRG thresholds (1/6, exact per Fedorenko) are NOT what the discrete BD simulations see:
  BD shows a splitting of alpha and alpha_s at theta_th = 0.25 +- 0.03 (Alés & López) and a change of
  height-distribution shape at theta_c ~= 0.2 (Jiang & Xia).
- Squizzato & Canet's NPRG finds SR-KPZ stable up to 1/6 in d = 1, in direct disagreement with
  Alés & López's numerics, and they say so explicitly (§2.4), attributing the discrepancy to
  finite-size/discretization effects.
- Alés & López conclude that no existing theory describes the faceted/anomalous-scaling regime, and that
  the DRG's infinitely many singularities signal its breakdown there.

---

## CONSOLIDATED NUMBERS

d = 1 unless stated. "computed" = my arithmetic from the quoted closed-form formula; "quoted" = printed
in a source (URLs in the sections above).

### A. Analytic branches (closed forms)

| Branch | theta range | chi = alpha | z | beta = chi/z |
|---|---|---|---|---|
| SR-KPZ (exact) | 0 <= theta <= 1/6 | 1/2 | 3/2 | 1/3 |
| Medina one-loop (Eq. 5) = Fedorenko Eqs. (79)-(80) = Katzav-Schwartz explicit SCE branch | theta > 1/6 | (1+4t)/(3+2t) | 2 - chi = 5/(3+2t) | (1+4t)/5 |
| Medina numerical fit (Eq. 7; "modified DRG") | theta >~ 1/6 | 1.69t + 0.22 | (2 chi+1)/(1+2t) | chi/z |
| Fedorenko two-loop | theta > 1/6 | - | 5/(3+2t+2 psi), psi <= 0, with 5/(3+2t) <= z <= 3/2 | - |
| Squizzato-Canet NLO NPRG | theta > 1/6 | 1.43t + 0.26 (fit) | printed: 1.49(4) @ 0.24, 1.42(5) @ 0.47; decreasing | not reported |
| Ma & Ma Flory/SA | all theta >= 0 | (1+2t)/(2+t) | 3/(2+t) | (1+2t)/3 |
| Katzav-Schwartz SCE, dominant (transcendental) branch | all theta >= 0 | no closed form (numerical, their Fig. 2) | no closed form; decreasing in t | - |
| Katzav-Schwartz cross-over between branches | - | - | - | phi_c = 1/6 (d = 1) |

### B. Values at selected theta (theory; computed unless marked quoted)

| theta | SR-KPZ (t <= 1/6) | Medina 1-loop chi/beta/z | Fedorenko z-window | Squizzato NPRG chi(fit)/z | Ma&Ma chi/beta/z |
|---|---|---|---|---|---|
| 1/6 = 0.1667 | 0.5 / 0.3333 / 1.5 | 0.5 / 0.3333 / 1.5 | [1.5, 1.5] | 0.4983 / 1.5 | 0.6154 / 0.4444 / 1.3846 |
| 0.20 | 0.5 / 0.3333 / 1.5 | 0.5294 / 0.36 / 1.4706 | [1.4706, 1.5] | 0.546 / - | 0.6364 / 0.4667 / 1.3636 |
| 0.24 | 0.5 / 0.3333 / 1.5 | 0.5632 / 0.392 / 1.4368 | [1.4368, 1.5] | 0.6032 / 1.49(4) (quoted) | 0.6696 / 0.5089 / 1.3393 |
| 0.25 | 0.5 / 0.3333 / 1.5 | 0.5714 / 0.40 / 1.4286 | [1.4286, 1.5] | 0.6175 / - | 0.6667 / 0.5 / 1.3333 |
| 0.30 | 0.5 / 0.3333 / 1.5 | 0.6111 / 0.44 / 1.3889 | [1.3889, 1.5] | 0.689 / - | 0.6957 / 0.5333 / 1.3043 |
| 0.35 | 0.5 / 0.3333 / 1.5 | 0.6486 / 0.48 / 1.3514 | [1.3514, 1.5] | 0.7605 / - | 0.7234 / 0.5667 / 1.2766 |
| 0.40 | 0.5 / 0.3333 / 1.5 | 0.6842 / 0.52 / 1.3158 | [1.3158, 1.5] | 0.832 / - | 0.75 / 0.60 / 1.25 |
| 0.47 | 0.5 / 0.3333 / 1.5 | 0.7107 / 0.552 / 1.2893 | [1.2893, 1.5] | 0.9321 / 1.42(5) (quoted) | 0.7696 / 0.6133 / 1.2411 |
| 0.50 | 0.5 / 0.3333 / 1.5 | 0.75 / 0.60 / 1.25 | [1.25, 1.5] | 0.975 / - | 0.80 / 0.6667 / 1.2 |

Squizzato & Canet's z is printed only at theta = 0.24 and 0.47; elsewhere it is readable only from their
Fig. 4(b) (hence "-"). Their paper never gives beta.

### C. Numerical (simulation) values, d = 1 and d = 2

| System | theta | alpha (global) | alpha_loc | alpha_s | z | beta |
|---|---|---|---|---|---|---|
| BD, d=1 (Alés & López) | 0.15 | 0.56(2) | - | 0.56(2) | 1.50 | - |
| BD, d=1 (Alés & López) | 0.47 | 1.02(1) | - | 1.23(2) | 1.50 | - |
| BD, d=1, columnar (Szendro et al. 2007) | 1/2 | 1.07(5) | - | 1.5(5) | - | - |
| BD threshold (Alés & López) | 0.25(3) | (splitting of alpha and alpha_s) | | | | |
| KPZ eq. threshold (Alés & López) | ~0.23 | (splitting of alpha and alpha_s) | | | | |
| BD & RSOS threshold (Jiang & Xia 2026) | ~0.2 | (change in height-distribution shape) | | | | |
| Continuum KPZ, d=1 (Hu, Hao, Xia 2023, PS scheme) | 0.30 | 0.73(3) | 0.788(5) | 0.828(6) | 1.25 (derived: alpha/beta) | 0.583(8) |
| Continuum KPZ, d=2 (Song & Xia 2016, quoted in arXiv:2306.06952) | 0.30 | 0.669 | - | - | 1.474 | - |
| Continuum KPZ, d=2, Gaussian noise (same source) | 0 | 0.385 | - | - | 1.571 | - |

### D. FIRM vs APPROXIMATE

FIRM (exact, or quoted directly from a source):
- theta_th = 1/6 in d = 1: exact crossover line 6 theta + 4 rho = 1 (Fedorenko Eq. 77, re-derived from his
  all-order-exact Eq. 74); same value from Medina et al. (as quoted) and from Squizzato & Canet's NPRG.
- SR-KPZ values chi = 1/2, z = 3/2 (beta = 1/3) for theta <= theta_th in d = 1.
- The exact LR relations z(1 + 2 theta) - 2 chi - d = 0 (S&C Eq. 6) / z_LR(1+2 theta) - 2 chi_LR + 2 rho = 1
  (Fedorenko Eq. 78).
- Fedorenko's bound 5/(3 + 2 theta) <= z_LR <= 3/2 for theta in [1/6, 1/2], and that it excludes Medina's
  "modified DRG" fit and the dominant second SCA solution.
- Squizzato & Canet: theta_th = 0.166 (d=1) and 0.346 (d=2); chi(theta) ~= 1.43 theta + 0.26 (NLO fit);
  z_LR = 1.49(4) at theta = 0.24 and 1.42(5) at theta = 0.47; chi_SR(d=2) ~= 0.375.
- Katzav & Schwartz: closed-form explicit branch z = (d+4)/(3+2 phi), Gamma = (d+4)(1+2 phi)/(3+2 phi),
  i.e. d = 1: alpha = (1+4 phi)/(3+2 phi), z = 5/(3+2 phi), beta = (1+4 phi)/5; Gamma_0(1) = 2,
  Gamma_0(2) = 2.59; phi_c(d=1) = 1/6.
- Alés & López and Jiang & Xia numbers in table C.

APPROXIMATE / EXTRAPOLATED (use with care):
- Medina's one-loop chi_LR = (1+4 theta)/(3+2 theta): "only valid for small theta close to the threshold".
- Medina's fit chi = 1.69 theta + 0.22: fit to numerically integrated truncated one-loop flow equations;
  ruled out by Fedorenko's stability bound.
- Squizzato & Canet's chi = 1.43 theta + 0.26 and z(theta): NLO derivative-expansion NPRG, fitted from
  figures, errors ~0.04-0.05 on z; their z values are consistent with their own chi fit through the exact
  LR relation (my check).
- Ma & Ma's Flory formulas: a scaling ansatz with no threshold; they violate the exact LR relation and
  Fedorenko's bound.
- Fedorenko's two-loop results: epsilon = 4 - d expansion extrapolated to epsilon = 3; his own words:
  "the large value epsilon = 3 of the expansion parameter describing the DP problem makes the estimation
  of psi highly unreliable".
- Everything from the discrete simulations: finite-size data-collapse estimates with the quoted errors;
  alpha_s is a genuinely new (non-standard) exponent, and its measured value depends on the fit window.
- The "no threshold" SCE scenario in d = 1: only numerically available (Fig. 2 of cond-mat/0309280), and
  the paper contains apparently contradictory statements about which of its two solutions dominates.

NOT OBTAINABLE from accessible sources:
- Any table of (theta, chi, beta, z) from Squizzato & Canet: no such table exists, and beta is never
  reported; in d = 2 the whole LR theta-dependence exists only as curves in their Fig. 6.
- Medina et al.'s paper itself (Phys. Rev. A 39, 3053, 1989): paywalled, not on arXiv; §1 is second-hand.
- Fedorenko's numerical z(theta) at epsilon = 3: figure-only (his Fig. 8), and he warns it is unreliable.
- The exact coefficient of psi in Fedorenko's Eq. (87) (text extraction gives "3 + 2 theta + psi";
  internal consistency requires "2 psi") — check the printed PRB page.
- Any source that treats theta < 0 explicitly (see §5; the closest is Katzav & Schwartz's sector-alpha
  statement for phi <= 0, which is about d > 2 / EW).
- Jiang & Xia's printed display equations (Eqs. 3-13) and any alpha_s values from that paper: the IOP
  full text renders them as images; the PDF is not fetchable with these tools.

---

## BibTeX

@article{Medina1989,
  author  = {Medina, E. and Hwa, T. and Kardar, M. and Zhang, Y.-C.},
  title   = {Burgers equation with correlated noise: Renormalization-group analysis and applications to directed polymers and interface growth},
  journal = {Phys. Rev. A}, volume = {39}, pages = {3053}, year = {1989},
  doi     = {10.1103/PhysRevA.39.3053}
}

@article{Squizzato2019,
  author        = {Squizzato, Davide and Canet, L{\'e}onie},
  title         = {Kardar-Parisi-Zhang equation with temporally correlated noise: A nonperturbative renormalization group approach},
  journal       = {Phys. Rev. E}, volume = {100}, pages = {062143}, year = {2019},
  doi           = {10.1103/PhysRevE.100.062143},
  eprint        = {1907.02256}, archivePrefix = {arXiv}, primaryClass = {cond-mat.stat-mech}
}

@article{Hanfei1993,
  author  = {Hanfei and Ma, Ben-Kun},
  title   = {Scaling analysis of Langevin-type equations},
  journal = {Phys. Rev. E}, volume = {47}, pages = {3738--3740}, year = {1993},
  doi     = {10.1103/PhysRevE.47.3738},
  note    = {cited in the literature as ``Ma \& Ma (1993)'' or ``Han \& Ma''}
}

@article{Katzav2004,
  author        = {Katzav, Eytan and Schwartz, Moshe},
  title         = {Kardar-Parisi-Zhang equation with temporally correlated noise: A self-consistent approach},
  journal       = {Phys. Rev. E}, volume = {70}, pages = {011601}, year = {2004},
  doi           = {10.1103/PhysRevE.70.011601},
  eprint        = {cond-mat/0309280}, archivePrefix = {arXiv}, primaryClass = {cond-mat.stat-mech}
}

@article{Fedorenko2008,
  author        = {Fedorenko, Andrei A.},
  title         = {Elastic systems with correlated disorder: Response to tilt and application to surface growth},
  journal       = {Phys. Rev. B}, volume = {77}, pages = {094203}, year = {2008},
  doi           = {10.1103/PhysRevB.77.094203},
  eprint        = {0712.0801}, archivePrefix = {arXiv}, primaryClass = {cond-mat.dis-nn}
}

@article{Ales2019,
  author        = {Al{\'e}s, Alejandro and L{\'o}pez, Juan M.},
  title         = {Faceted patterns and anomalous surface roughening driven by long-range temporally correlated noise},
  journal       = {Phys. Rev. E}, volume = {99}, pages = {062139}, year = {2019},
  doi           = {10.1103/PhysRevE.99.062139},
  eprint        = {1902.06674}, archivePrefix = {arXiv}, primaryClass = {cond-mat.stat-mech}
}

@article{Lam1992,
  author  = {Lam, Chi-Hang and Sander, Leonard M. and Wolf, Dietrich E.},
  title   = {Surface growth with temporally correlated noise},
  journal = {Phys. Rev. A}, volume = {46}, pages = {R6128}, year = {1992},
  doi     = {10.1103/PhysRevA.46.R6128}
}

@article{Szendro2007,
  author  = {Szendro, I. G. and L{\'o}pez, J. M. and Rodr{\'i}guez, M. A.},
  title   = {Localization in disordered media, anomalous roughening, and coarsening dynamics of faceted surfaces},
  journal = {Phys. Rev. E}, volume = {76}, pages = {011603}, year = {2007},
  doi     = {10.1103/PhysRevE.76.011603}
}

@article{Jiang2026,
  author  = {Jiang, Pei-Yan and Xia, Hui},
  title   = {Height distributions of discrete growth models in the presence of long-range temporal correlations},
  journal = {Communications in Theoretical Physics}, volume = {78}, pages = {105602}, year = {2026},
  doi     = {10.1088/1572-9494/ae8f20}
}

@article{Hu2023,
  author        = {Hu, Xiongpeng and Hao, Dapeng and Xia, Hui},
  title         = {Improved finite-difference and pseudospectral schemes for the Kardar-Parisi-Zhang equation with long-range temporal correlations},
  journal       = {Physica A}, volume = {619}, pages = {128744}, year = {2023},
  doi           = {10.1016/j.physa.2023.128744},
  eprint        = {2105.13528}, archivePrefix = {arXiv}
}

@article{ChangXia2022,
  author        = {Chang, Zhichao and Xia, Hui},
  title         = {Probability distributions for kinetic roughening in the Kardar-Parisi-Zhang growth with long-range temporal and spatial correlations},
  year          = {2022}, eprint = {2205.08483}, archivePrefix = {arXiv}, primaryClass = {cond-mat.stat-mech}
}

@article{SongXia2016,
  author  = {Song, Tian-Shu and Xia, Hui},
  title   = {Long-range temporal correlations in the Kardar-Parisi-Zhang growth: numerical simulations},
  journal = {J. Stat. Mech. Theory Exp.}, volume = {2016}, pages = {113206}, year = {2016},
  doi     = {10.1088/1742-5468/2016/11/113206}
}

@article{Strack2015,
  author        = {Strack, Philipp},
  title         = {Dynamic criticality far from equilibrium: One-loop flow of Burgers-Kardar-Parisi-Zhang systems with broken Galilean invariance},
  journal       = {Phys. Rev. E}, volume = {91}, pages = {032131}, year = {2015},
  doi           = {10.1103/PhysRevE.91.032131},
  eprint        = {1408.1405}, archivePrefix = {arXiv}
}

@article{Canet2025,
  author        = {Canet, L{\'e}onie},
  title         = {The non-perturbative sides of the Kardar-Parisi-Zhang equation},
  year          = {2025}, eprint = {2509.01472}, archivePrefix = {arXiv}, primaryClass = {cond-mat.stat-mech}
}
