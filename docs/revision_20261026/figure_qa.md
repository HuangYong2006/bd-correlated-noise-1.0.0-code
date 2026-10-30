# Revised figure contract and QA

Date 2026-10-26. Backend is Python/matplotlib, confirmed by figure_backend.json. All four figures have been rendered from the completed revision data and visually inspected panel by panel. Only the previous 7 pt typography, restrained colours and 178 mm width were inherited. No old curve, threshold shading or old physical conclusion was reused.

## Sources and reproduction

- Plotting source is code/revision_figures.py. The existing code/figures_ctp.py keeps its historical helper functions, while its final main block calls the new entry point.
- Input summary is data/analysis/revision_20261026/summary.json, SHA256 `57ba7392911c2677fb50d3a9db82da3b1d9b593a711ddc10bd60598622ba6091`.
- Noise input is data/analysis/revision_20261026/noise_validation.json, SHA256 `69abe021091762baf559b92b21a0c9a211fa94003b672e30be803ac22bc3b8fe`.
- Rendered outputs are figs/revision_20261026/fig1_generator, fig2_exponents, fig3_fixed_theta and fig4_morphology, each in PDF, SVG and PNG. figure_sources.json records source hashes, the full protocol, raw files, representative-profile indices and output paths.
- The declared simulation tags are rev20261026_shared, rev20261026_powerlaw and rev20261026_fgn. There are 64 independent runs per cell, five sizes from 32 through 512 and seven theta values. The script checks the raw morphology files belong to those exact campaigns and that the spectral files carry noise_layout="spectral-run-time-site-v2".

```powershell
& 'D:/TOOLS/anaconda3/python.exe' code/revision_figures.py
# The historical CLI also routes to exactly the new renderer.
& 'D:/TOOLS/anaconda3/python.exe' code/figures_ctp.py --summary data/analysis/revision_20261026/summary.json --outdir figs/revision_20261026
```

Both commands were executed successfully with the real inputs. No placeholder or demonstration data were generated.

## Evidence and panel contract

| Panel | Question and displayed data | Centre, uncertainty and units | Visual QA |
|---|---|---|---|
| 1a | Does production fGn at theta=0.25, L=64 follow its finite-record temporal target? | Eight independent runs, binary correlation mean ± run SE. Dashed target is the finite-period, DC-removed production spectral covariance passed through Gaussian arcsine, not the infinite-fGn covariance. | Target and production nearly overlap as expected; markers and target remain distinguishable, legend clears data. |
| 1b | Does the old layout introduce spatial correlation? | Current and legacy layouts use the same underlying theta=0.25, L=32 fGn samples. Eight-run correlation means ± SE; zero reference. | Matching sizes are stated, labels clear the curves. Small current errors are visible near zero at sufficient zoom. |
| 1c | Do shared and independent AR innovations yield the same theta=0 covariance? | Exact stationary normalised raw-field covariances, N=48. No binary conversion, no uncertainty for these deterministic sums. This is not the simulator's explicit iid theta=0 branch. | Distinct solid and dashed lines, positive log inputs, no theoretical phase boundary. |
| 2a | How do finite-size late-time width slopes vary with theta in the three new campaigns? | All 21 generator/theta fits; all five sizes and all 64 runs per cell. Error bars are supplied run-bootstrap SD from 500 resamples. Only alpha=1/2 is a reference. | All three trajectories rise; no old flat-KPZ narrative, threshold region or extrapolated theory curve remains. |
| 2b | How do adjacent-size estimates vary at theta=0.25? | All four size pairs for each generator, supplied run-bootstrap SD. The x coordinate is the summary's geometric pair midpoint. | Every pair and error bar is shown, no clipping, shared legend applies to both panels. |
| 3a | How does the local flat-pair fraction change at L=256? | All seven theta cells for all three generators, all 64 runs, supplied run-bootstrap SD from 500 resamples. | Curves and legends are separated. The figure makes no assertion about universality from this local statistic. |
| 3b | How does the local-pmf distance vary with size for fGn? | All five sizes at theta=0.1, 0.25 and 0.4, all 64 runs. dH_err is the supplied outer independent-run bootstrap SD with 200 draws; conditional null scatter is not substituted. | All three series and error bars are visible, including the low signal; no crossing or significance marker. |
| 4a | What do representative final profiles look like? | One illustrative final profile per generator at theta=0.4, L=256, selected by final width nearest the 64-run median. Mean-centred, then vertically translated without rescaling. Dotted baselines and a 10-height-unit bar show the shifts and physical scale. | Full width is shown, profiles and scale bar do not collide. This panel is explicitly illustrative. |
| 4b | What is the ensemble spatial spectrum? | All 64 final profiles for each generator. S(k)=abs(rFFT(h-mean(h)))²/L, all positive modes. Mean ± one SE, using run SD/sqrt(64). | All three spectra and uncertainty bands are shown with no smoothing, clipping or artificial positive floor. |

Figure archetypes are quantitative grids, with an illustrative-profile plus quantitative-spectrum composition for Fig. 4. A rise in a local statistic is not labelled as a necessary or sufficient universality-crossover test. Finite-size slopes are not labelled asymptotic exponents.

## Data inclusion and transformations

No simulation run is excluded from a quantitative panel. Fig. 2 uses all 105 cell inputs through the supplied fits. Fig. 3 intentionally selects L=256 for its theta scan and three predefined theta values for its size scan. Fig. 1 intentionally shows the prescribed theta=0.25 covariance validation and theta=0 AR audit; the JSON retains all eight production-validation conditions and per-run values.

The zero-based representative-profile indices are 23 for shared AR, 5 for power law and 1 for fGn. These choices are reproducible nearest-median-width selections, not convenient examples chosen after inspection. Their raw files are respectively

- data/raw/rev20261026_shared/bd_p0d4000_L256.npz
- data/raw/rev20261026_powerlaw/bd_p0d4000_L256.npz
- data/raw/rev20261026_fgn/bd_p0d4000_L256.npz

Fig. 4 drops only the zero Fourier mode after subtracting each profile's mean. No positive-frequency bin is omitted or pooled. All quantities passed to logarithmic y axes are checked to be positive, including the lower spectrum uncertainty limit. Alpha and local statistics are taken directly from the new summary rather than recalculated with another convention in the plotting code.

## Export and verification

PDF page dimensions were read back independently with pypdf. All widths are 178.000 mm. Figure heights are 80, 78, 78 and 84 mm, respectively. PNGs are 600 dpi, with dimensions 4204×1889, 4204×1842, 4204×1842 and 4204×1984 pixels. SVG text remains text; PDF uses editable embedded TrueType fonts.

The nature-figure static preflight reported 17 PASS, 3 WARN and 0 FAIL. Warnings were inspected rather than ignored

- The static detector cannot evaluate the expression 178/25.4 and incorrectly reports 4521.2 mm. Actual PDF widths were measured as 178.000 mm.
- The requested raster deliverable is PNG; TIFF was not requested. PDF and SVG are the vector deliverables.
- The conservative mathtext warning combines the smallest legend size with a larger axis math label. Actual PDF font audits passed for every text run. Minimum glyph sizes are 5.04 pt, 5.25 pt, 5.25 pt and 5.04 pt for Figs. 1 through 4.

Every panel and each assembled PNG was inspected. Axis labels, titles, panel letters and legends are within the page with no clipping or collisions; markers and uncertainty bounds are visible. Direct and bridged rendering preserve the same 178 mm canvas because the new style explicitly resets the historical tight-bbox setting.

The full repository suite was run after plotting and passed `48 passed in 2.50s`. `git diff --check` passed for the plotting files. This QA verifies data routing, plotting and rendering; it does not establish asymptotic universality, a critical theta, or the implementation used by external authors.

## Final root checks

The root agent inspected all four regenerated PNGs and confirmed the panel
meanings and uncertainty definitions against the final captions. Figure 2's
ordinate was clarified to finite-size slope alpha_fit. SVG export now strips
trailing whitespace; this does not alter the paths or rendered appearance.

The static figure preflight reports 17 PASS, 3 WARN, 0 FAIL. Its width parser
misreads the 178/25.4 expression; actual PDF widths are 178 mm. PNG previews plus
vector PDF/SVG satisfy the present PRE source package, so TIFF is not needed.
Mathtext shrinkage was checked on the PDFs: the minimum text runs are 5.04 pt
in Figures 1/4 and 5.25 pt in Figures 2/3, with zero runs below 5 pt. These
resolve the three warnings without changing the scientific content.

## Post-review corrections

The two nonblocking manuscript-audit suggestions were resolved without changing
quantitative inputs. The figure labels now consistently say fGn spectrum, and
the separate theory note states the positive-integer lag domain of its fGn
covariance expression. The figures and Overleaf package were rebuilt afterward.
The frozen final review reports are not rewritten to hide the original findings.
