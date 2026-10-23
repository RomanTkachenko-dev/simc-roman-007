Combined signal update with corrected rho weights, 2026-10-05

Main output: combined_signal_weighted_rho.png / .svg
All seven templates use the original common acceptance cuts, 100 bins
over missing mass 0-2 GeV. The black Pythia curve and its display multiplier
are unchanged. Only the rho input template was replaced: every accepted
rho event now contributes its corrected positive Weight instead of 1.
All seven nonnegative amplitudes were then refitted with fixed shapes.
No smoothing, missing-mass shifts, or extra background terms were added.
The prior raw templates and prior fit were independently reproduced first.

rho input: ../rho_jacobian_fix_20261005/after/rho_Roman_007.root
573 accepted events; sum Weight=0.0001995499496037067.
Weighted templates: omega (still provisional), rho, 2pi, 3pi.
Count templates: phi, eta, eta-prime.

Solver: enumeration of all 128 possible active parameter sets, solving
linear least squares for each feasible set and selecting the global minimum.
The objective is identical to the old NNLS/ROOT fit: observed Pythia count
variance, with variance 1 in empty bins. Verified by old-result reproduction,
KKT optimality and synthetic closure. No new ROOT Minuit fit was performed.
results.root stores input and fitted histograms, the sum, Pythia, residuals,
and JSON metadata. It does not contain a ROOT minimizer object/canvas.

Fit diagnostics:
chi2 changed from 874426.1773 to 870380.6034, with nominal ndf=93.
Large residuals remain. This fit is not statistically adequate under its
Pythia-only error model, which neglects finite MC shape uncertainty.
Omega count equivalent changed from 1173838.82 to 1163130.62 (about -0.91%).
Rho count equivalent changed from 1024908.14 to 1081758.24 (about +5.55%).
These numbers are diagnostic Pythia count equivalents, not physical yields
or cross sections. No absolute luminosity normalization has been inferred.

rho_template_shape_comparison.png compares unit-area rho shapes before/after
weighting. Its normalization is for shape comparison only.

Omega status:
The user is unsure whether the intended basis is an updated SIMC model or
Henry's version. No replacement omega function/formula/revision is supplied.
The existing omega sigcm is negative for all 595 accepted baseline events.
The corrected rho Jacobian does not alter the omega model. Its old weighted
template is retained and explicitly labelled provisional in the combined plot.

Henry reference:
The saved Henry SIMC examples here concern phi and J/psi, not a validated
omega sample. The supplied /path/to/EIC/Reference/plot_exp_1over35.py assumes
dsigma/d|t|=A exp(-5|t|), with A=10,20,30,40 for W intervals
1.80-2.00, 2.00-2.05, 2.05-2.10, 2.10-2.40 GeV, respectively.
It does not extract these central cross-section values from its fitted yields.
henry_omega_reference.png / .svg reproduces only those assumed curves.
henry_omega_reference.json includes the values, exact source hash and limits.
No statistical uncertainties or absolute SIMC cross-section claim are added.

Next required physics input:
Bill needs to identify the omega parameterization valid for the chosen
W,Q2,t/u,epsilon and the intended Jacobian. Henry's curves can be used as a
clearly stated reference assumption, but not silently adopted as a validated
SIMC model. Once specified, rerun omega and recompute the physical prediction.
