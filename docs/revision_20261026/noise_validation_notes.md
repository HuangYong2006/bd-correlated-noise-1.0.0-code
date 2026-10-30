# Production-noise validation notes

The executable check is code/validate_production_noise.py. It validates the final
corrected production interface, not a separate one-dimensional generator. The
saved data are data/analysis/revision_20261026/noise_validation.json.

The eight declared cells are powerlaw/fgn times theta=.25/.40 times L=32/64,
with T=65536 and eight independent runs. Temporal lags are 1,2,4,8,16,32,64,128;
spatial lags are 1,2,4,8. Every cell passes the declared engineering check of
five run standard errors plus absolute tolerance .003. Individual-run values
are preserved, so the SE definition can be checked directly.

The target is the actual finite-period Gaussian covariance obtained by inverse
FFT of the production spectrum, with DC removed and the existing alias cutoff.
Gaussian sign covariance follows the arcsine identity. The target is not the
infinite-record fGn covariance. This boundary matters especially at strong
correlation. The script records both conventions for reference.

The legacy layout control uses exactly the same fGn theta=.25,L=32 samples,
transposed to (run,site,time) and reshaped to the original (run,time,site) shape.
It fails both time and space checks. The measured correct lag-one temporal
correlation is .270413 against finite-period target .270190; correct spatial
lag-one correlation is .000012. The faulty mapping yields temporal .039991 and
spatial .263201. This is a controlled test of this repository's layout fault,
not a statement about another author's code.

The ar_zero fields are deterministic stationary raw covariances at theta=0,
N=48, a=6, B=2. They retain cross terms for common innovations and are compared
with independent innovations. They are not binary-field covariances, and the
Gaussian arcsine identity is not imposed on the uniform-driven AR input.
They are also not the production simulator's explicit iid shortcut at theta=0.

The executable stores seeds, per-run arrays, definitions, tolerances and source
hashes. These checks cover the declared cells only; they do not prove every
possible record length or an asymptotic universality class.
