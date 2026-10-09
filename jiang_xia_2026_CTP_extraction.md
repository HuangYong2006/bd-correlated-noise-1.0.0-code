# Extraction report — Jiang & Xia, Commun. Theor. Phys. 78, 105602 (2026)

DOI: 10.1088/1572-9494/ae8f20 · "Height distributions of discrete growth models in the presence of
long-range temporal correlations" · Peiyan Jiang, Hui Xia (China Univ. of Mining and Technology)

## Sources actually used (all fetchable)

| Purpose | URL | Result |
|---|---|---|
| Full text (HTML, accepted manuscript) | https://beta.iopscience.iop.org/article/10.1088/1572-9494/ae8f20/meta | HTTP 200, full text retrieved |
| Same via W3C html-to-text (full text in one blob) | https://www.w3.org/services/html2txt?url=https%3A%2F%2Fbeta.iopscience.iop.org%2Farticle%2F10.1088%2F1572-9494%2Fae8f20%2Fmeta | HTTP 200, full text |
| Publisher-of-record landing page | https://iopscience.iop.org/article/10.1088/1572-9494/ae8f20 | Bot-gated (perfdrive redirect via web_fetch); via html2txt returns only title+abstract+solid "Access this article" paywall block |
| Complete reference list (64 refs, with DOIs) | https://api.crossref.org/works/10.1088/1572-9494/ae8f20 | HTTP 200 JSON |
| OA/licence + content locations | https://api.openalex.org/works/doi:10.1088/1572-9494/ae8f20 | CC-BY-NC-ND, hybrid OA, publishedVersion |
| Semantic Scholar record | https://api.semanticscholar.org/graph/v1/paper/DOI:10.1088/1572-9494/ae8f20?fields=title,externalIds,openAccessPdf,abstract | CCBYNCND, OA PDF = IOP pdf only |
| arXiv check | https://export.arxiv.org/api/query?search_query=au:%22Hui+Xia%22+AND+abs:%22temporal%22&max_results=40 | **No arXiv version of this paper exists** (only earlier Xia-group papers) |

Raw spill file with the complete extracted article text (greppable):
`C:\Users\MSI\AppData\Local\Temp\dsh-spill-Ot5J1V\session-65e93bcf52a4\558bf999c28b-web_fetch.txt`

## 1. Models and noise generation — what is verifiable

**Noise correlator** (§2, printed inline as text, verbatim):
> "one needs to generate a long-range correlated sequence that satisfies the scaling relation:
> $C(\tau) = \langle \eta_{\theta}(t)\eta_{\theta}(t+\tau) \rangle \sim \tau^{2\theta-1}$"

**Generator**: "we adopt the fast fractional Gaussian noise (FFGN) technique [61] and the binary mapping
method [62] for this purpose", with the weight function (Eq. 6) parameters given in text as
"$B = 2, a = 6, u_{n} = aB^{-n}, r_{n} = \mathrm{e}^{-u_{n}}$", number of components $N = 40$,
and FFGN becoming numerically unstable for $\theta > 0.4$.

**BOUNDARY CONDITIONS** (verbatim): "starting from an initially flat interface, under periodic boundary
conditions where the leftmost and rightmost columns are treated as nearest neighbors." Substrate size $L = 1024$.

**NOT OBTAINABLE**: the printed formulas of Equations (3), (4), (5), (6), (7), (8), (9), (10), (11).
In the HTML these display equations are rendered client-side as images; the extractor sees only the
literal placeholder text `Equation (10)`, etc. The PDF is unreadable by web_fetch ("unsupported content
type"), pwsh has no network, and OpenAlex's GROBID/PDF content endpoints now require an API key
(`https://content.openalex.org/works/W7170175408.grobid-xml` → HTTP 401 "API key required").
Prose descriptions that DO survive (verbatim):
- BD: "all incident particles are released from a randomly selected position far above the aggregate,
  follow a vertical trajectory until they reach the aggregate front, and stick upon the first contacting
  nearest neighbor occupied site"
- RSOS: "each incident particle is incorporated at the top of the target column, provided that height
  differences between the selected and neighboring columns satisfy the condition $|\Delta h| \leq 1$.
  When this critical height difference is exceeded, the deposition attempt is rejected."
- The temporally correlated BD model is cited to refs [48],[60],[62]; the RSOS rule is "given by the
  following mathematical expression" (Eq. 11, not extractable).
- Sign convention: "BD and RSOS models correspond to positive and negative values of the nonlinear
  coefficient $\lambda$"; RSOS results are shown after the transformation $h \to -h$.

## 2. Quantitative findings vs theta — NO TABLES EXIST

The article contains **no tables at all** (only a nav link labelled "Tables"). All quantitative results
are in figures; the only numbers printed in text/captions are:

- $\theta \in [0, 0.4]$ (simulation range); stable/simulations reliable within this range.
- $\theta_\mathrm{c} \approx 0.2$ in the growth regime (text and abstract): "for $\theta < \theta_\mathrm{c}$
  the system exhibits a stable plateau consistent with the TW-GOE distribution"; for $\theta > \theta_\mathrm{c}$
  "a crossover toward a new universality class".
- BD: height distributions "qualitatively similar to the TW-GOE form" for $0 \leq \theta \leq 0.2$;
  "for $\theta \geq 0.25$, the right tail of distributions exhibits a progressive deviation"; $S$ and $K$
  "remain stable and approach the exact values for TW-GOE up to $\theta = 0.2$ ... when $\theta > 0.2$,
  both $S$ and $K$ exhibit a pronounced increase".
- RSOS: "$S$ and $K$ ... As $\theta$ increases to $0.2$, both estimated values exhibit only a slight
  increase and remain close to those of TW-GOE. Notably, a pronounced increase is observed in the range
  $0.2 \leq \theta \leq 0.4$."
- Reference lines in captions: TW-GOE exact values $S = 0.2934$, $K = 3.1652$ (figs 3, 5); BR-$F_0$ exact
  values $S = 0.3594$, $K = 3.2891$ (figs 9, 11); Gaussian $S = 0$, $K = 3$.
- Saturation-regime height-increment crossover: BD — $S,K$ decrease for $\theta \in [0,0.3]$, "pronounced
  increase in the regime $0.3 < \theta < 0.4$"; RSOS — decrease for $0 < \theta < 0.25$, sharp increase
  for $\theta > 0.25$.
- Figure sampling: $\theta = 0.0, 0.1, 0.2, 0.25, 0.3, 0.4$ (figs 2, 4, 8, 10); $\theta = 0.05, 0.2, 0.4$
  (fig 13). Morphology snapshots at $t = 200$, $t = 2000$, $t = 2^{18}$.
- Growth exponent: "the growth exponent $\beta$ is found to increase with $\theta$, in qualitative
  agreement with previous studies on temporally correlated BD and KPZ growth" — **no $\beta$ values**;
  $\beta(\theta)$ exists only as an inset in figure 1.

**Exponents $\alpha$, $\alpha_{loc}$, $\alpha_s$, $z$: NOT reported anywhere in this paper.** The only
exponent values present are the uncorrelated KPZ reference values in the introduction ($\alpha = 1/2$,
$\beta = 1/3$ for (1+1)-dimensional KPZ) and the definitions $z = \alpha/\beta$ and the identity
$\alpha + z = 2$.

## 3. Comparison with continuum KPZ / theory (verbatim, all from the beta full text)

- Medina et al.: "Specifically, Medina et al proposed that when the correlation exponents are below a
  critical threshold $\theta_\mathrm c = 0.167$, the growth exponents remain independent of the
  correlation exponents through dynamic renormalization group [41]." Ref [41] = Medina et al.,
  Phys. Rev. A **39**, 3053 (1989), DOI 10.1103/PhysRevA.39.3053. **No numbers attributed beyond 0.167.**
- Han & Ma: "In contrast, Han and Ma argued that the critical exponents vary linearly with the correlation
  exponent by the Flory-like scaling approach [11]." Ref [11] = Hanfei, "Scaling analysis of Langevin-type
  equations", Phys. Rev. E **47**, 3738 (1993), DOI 10.1103/PhysRevE.47.3738. **No numbers.**
- Katzav & Schwartz: "Through self-consistent expansion, Katzav and Schwartz found a continuously varying
  form for the exponents [42]." Ref [42] = Katzav & Schwartz, Phys. Rev. E **70**, 011601 (2004),
  DOI 10.1103/PhysRevE.70.011601. **No numbers.**
- Squizzato & Canet (PRE **100**, 062143 (2019), DOI 10.1103/PhysRevE.100.062143) appears as reference
  **[45] in the reference list but is never cited individually in the body**; the only plausible mention
  is the range citation: "The validity of these theoretical predictions has been qualitatively verified by
  numerical simulations of the KPZ equation with long-range correlations [43–47]." **No numbers at all.**
- Fedorenko, Phys. Rev. B **77**, 094203 (2008): **not cited anywhere.** It does not appear in the
  64-item reference list (checked in the Crossref reference array).
- Comparison with continuum KPZ with correlated noise: "Our results are further compared with recent
  findings on height fluctuations in the continuum KPZ equation driven by long-range temporally correlated
  noise [58], as well as with earlier studies of macroscopic dynamic behavior [46, 59, 60]."
  Ref [58] = Chang & Xia, Commun. Theor. Phys. **78** (2026), DOI 10.1088/1572-9494/adfa78;
  [59] = Rodríguez-Fernández et al., Phys. Rev. E **110**, 024104 (2024);
  [60] = Wang et al., Phys. Rev. E **111**, 024124 (2025); [46] = Alés et al., Phys. Rev. E **99**, 062139 (2019).
- Threshold comparison (conclusion, verbatim): "We find a critical threshold $\theta_\mathrm c = 0.2$ in
  discrete growth models, consistent with observations for the continuous KPZ equation [58, 59]. However,
  this value differs from the threshold $\theta = 0.25$ derived from macroscopic dynamic phenomena
  [46, 60], suggesting that height fluctuations are more sensitive to long-range temporal correlations."
- Statement of absence of theory: "At present, a precise theoretical prediction for the KPZ universality
  class with long-range temporal correlations is still lacking."

## 4. Anomalous scaling / faceting

- The words $\alpha_{loc}$ and $\alpha_s$ and "spectral roughness exponent" **never appear** in this
  paper; there is no claim about $\alpha_s > 1$.
- Anomalous scaling is only referenced, not measured: "studies focusing on dynamic properties, including
  anomalous scaling [46, 59] and emergent phenomena [60], have reported a critical value $\theta = 0.25$."
- Galilean-invariance breaking: "this coupling breaks Galilean invariance within the KPZ framework,
  manifesting as a violation of the scaling identity $\alpha + z = 2$ [46–48, 59]."
- Faceting: "The asymptotic noise correlator of temporally correlated KPZ class in the limit
  $\theta \to 0.5$ exhibits theoretical consistency with columnar-disordered KPZ class [63], and their
  interface morphologies are qualitatively similar in a faceted pattern [59]." Morphology: for
  $\theta = 0.4$ "distinct local faceted morphology emerges" at $t = 2000$ and "the interface ultimately
  develops global faceted patterns" at $t = 2^{18}$; $\theta = 0.05, 0.2$ remain self-affine.
- Conclusion caveat: "due to the limitations of the FFGN method, we do not obtain sufficiently conclusive
  results in the limit $\theta \to 0.5$."
