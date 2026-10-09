# Cover letter — Physical Review E

<!--
HOW TO USE THIS FILE
  1. Fill in the four bracketed items: date, corresponding e-mail, and the
     referee section at the end.
  2. Paste the body into the APS submission server's cover-letter box, or
     upload it as a separate PDF.
  3. Remember that the Data Availability Statement has to be entered in the
     submission system as well as appearing in the paper.
-->

[Date]

To the Editors of *Physical Review E*:

Dear Editors,

We submit for your consideration the manuscript **"Noise-generator bias in
temporally correlated ballistic-deposition simulations"** for publication as a
Regular Article.

## Context of the results

Numerical studies of kinetic roughening driven by noise with long-range temporal
correlations disagree about the exponent theta_c at which Kardar-Parisi-Zhang
scaling gives way to a long-range regime; the reported values span from
theta_c -> 0 to theta_c = 0.25 +/- 0.03. All of these studies generate their noise
by thresholding a fast fractional Gaussian noise (FFGN) sequence, so the nominal
correlation exponent of the generator is identified with the physical exponent of
the noise correlator. Our manuscript examines whether that identification holds.

## Summary of key findings

1. **The FFGN construction documented in the studies we could verify does not
   approach white noise continuously as the nominal exponent goes to zero.** Its
   normalised correlation amplitude is essentially independent of theta, whereas
   the exact fractional-Gaussian-noise amplitude 2*theta*(theta + 1/2) vanishes;
   at nominal theta = 0 the generated noise therefore retains a C(tau) ~ tau^-1
   tail instead of being white. We verify by converging the generator in the
   number of components that this is intrinsic to the construction and not a
   truncation artefact.

2. **A shape-matched calibration performed on the binarised field** — which is
   the calibration the discrete model requires — **maps the entire published
   theta axis onto effective values at or above 0.19**, i.e. above the
   renormalisation-group threshold 1/6. The family as documented never samples the
   short-range regime.

3. **Inside the model this changes the physics that is measured.** At nominal
   theta = 0 it produces a 23% discontinuity of the structural order parameter,
   inverts the theta-dependence of the interface morphology, and inflates the
   apparent variation of the roughness exponent by a factor of eight.

4. **A theta-continuous generator that carries the correct correlation
   amplitude does not reproduce the reported onset behaviour** within the system
   sizes we can reach, and the deviations it shows behave as finite-size
   corrections that shrink with L rather than as asymptotic behaviour.

5. **The standard construction for locating a crossover is confounded.** A
   fixed-significance crossing of a distributional distance necessarily drifts
   when the detection floor depends on system size; we measure that drift,
   identify its origin, and show that testing scale dependence at *fixed*
   parameter is the appropriate discriminator. We report the result of that test
   and a feature-ablation analysis of the two diagnostics used.

## What we do not claim

This matters for how the paper should be read, so we state it explicitly here as
well as in the manuscript. **We do not report a measured theta_c.** Within the
accessible window L = 32 to 512, no crossover is resolved for theta <= 0.25 with
a generator of correct correlation amplitude, and the renormalisation-group
prediction 1/6 is neither confirmed nor contradicted by our data. Our claim is
narrower: the reported thresholds of 0.20 and 0.25 can be accounted for by the
generator bias within the accessible size range, and the onset construction
cannot locate a crossover in the presence of a system-size-dependent detection
floor. We also do not claim that all published work uses an identical
implementation; Table I states, reference by reference, which implementation
details we could establish from the published descriptions and which we could
not.

## Scope

The study is confined to the 1+1-dimensional parallel-update ballistic-deposition
model, to system sizes L = 32 to 512, and to generator implementations whose
published descriptions we were able to verify. We make no claim about
implementations whose parameters are not documented, and no claim about other
universality classes.

## Why we believe this is appropriate for Physical Review E

The result is methodological and numerical rather than a new measurement: it
identifies a reproducible bias in a widely used noise construction, and it shows
that a standard finite-size diagnostic can produce a spurious crossover when the
noise generator is not validated. We think it is of interest to the readership
that works on kinetic roughening, on finite-size scaling methodology, and more
generally on the reliability of numerical evidence for universality-class
changes. The paper combines exact renormalisation-group relations, an explicit
generator construction and calibration, two complementary distributional
diagnostics with null models, and a fixed-parameter finite-size test; the
negative and methodological parts are stated as such rather than presented as a
determination of the transition.

## Statements

- A **Data Availability Statement** is included in the manuscript: the simulation
  data, analysis scripts and simulator are available from the corresponding
  author on reasonable request.
- The manuscript is **not** part of a joint submission, and it has **not** been
  submitted previously to any Physical Review journal.

## Referees

<!-- Fill in, or delete this section if you have no preference. -->

Recommended referees:

- [name, affiliation, e-mail]
- [name, affiliation, e-mail]
- [name, affiliation, e-mail]

Referees we would prefer to be excluded:

- [name, affiliation, and a one-line reason if you wish]

Thank you for considering our manuscript.

Sincerely,

Yong Huang
College of Physics, Chengdu University of Technology, Chengdu, China
[corresponding e-mail]
