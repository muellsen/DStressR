# Future Extension: Dynamic Promoter-Activity Response Modeling

The PFIbox/French data suggest a natural extension of DStressR from scalar
reporter-growth summaries to dynamic promoter-activity response construction.
The present package model starts from positive scalar summaries and defines

\[
r_i=\log L_i-\alpha_{a_i}\log G_i,
\]

where \(L_i\) is the reporter summary, \(G_i\) is the matched growth or
biomass summary, and \(\alpha_{a_i}\) is a reporter-specific response exponent.
This is appropriate for Lux/OD or fluorescence/OD AUC summaries, but it does
not fully represent the established Zaslaver/French promoter-activity workflow.

For fluorescence time series, the original promoter-activity construction is
closer to a curve-level normalization. Let \(L_i(t)\) denote reporter
fluorescence, \(B_i(t)\) promoterless/background fluorescence, and \(G_i(t)\)
a biomass or constitutive-reporter proxy. The original workflow corresponds to

\[
q_i(t)=\frac{L_i(t)-B_i(t)}{G_i(t)}
\]

followed by

\[
PA_i(t)=\frac{d}{dt}q_i(t).
\]

Thus the fixed-exponent workflow is the special case \(\alpha_{a_i}=1\) of

\[
q_i^{(\alpha)}(t)=
\frac{L_i(t)-B_i(t)}{G_i(t)^{\alpha_{a_i}}},
\qquad
PA_i^{(\alpha)}(t)=\frac{d}{dt}q_i^{(\alpha)}(t).
\]

A future DStressR extension could estimate \(\alpha_a\) from control traces
rather than from scalar AUCs. One possible model is

\[
\log\{L_i(t)-B_i(t)+\epsilon\}
=
\eta_{a_i}
+
\alpha_{a_i}\log\{G_i(t)+\epsilon\}
+
h_{a_i}(t)
+
e_i(t),
\]

where \(h_a(t)\) captures the average reporter-specific time profile in
controls or is restricted to a stable time window. After estimating
\(\alpha_a\), DStressR could construct

\[
y_i^{PA,\alpha}
=
\frac{1}{|\mathcal T|}
\int_{\mathcal T} PA_i^{(\alpha)}(t)\,dt
\]

and pass this scalar response into the usual DStressR effect-estimation and
testing layer.

This extension is attractive because it connects DStressR directly to the
established promoter-activity tradition while retaining the current downstream
inference framework. The main statistical choices are the smoothing method, the
time window \(\mathcal T\), weighting over time, background-trace calibration,
quality filtering for low reporter/background-corrected signal, and diagnostics
for stability of \(\widehat\alpha_a\). The original French/PFIbox workflow is
therefore not contradicted by the DStressR response model; rather, it points to
a dynamic response-construction module in which fixed \(\alpha_a=1\) is the
published baseline and estimated \(\alpha_a\) is a testable extension.

## Future Extension: Faster Interaction Models

The current full interaction path in `fit_destress(..., interaction = TRUE)`
uses ordinary dense linear models. This is statistically transparent and useful
for small to moderately sized interaction problems, but it is not the right
computational backend for very large screens. The model matrix contains
reporter main effects, perturbation main effects, technical covariates, and
reporter-by-perturbation interaction terms. Its interaction block scales as
approximately

\[
n(p-1),
\]

where \(n\) is the number of reporters and \(p-1\) is the number of
non-reference perturbations. This is feasible for tall data with few
perturbations, such as the drug-level PFIbox analysis, but it becomes expensive
for wide compound screens.

Several implementation routes should be considered before exposing the full
interaction model as a routine large-screen default:

- Use sparse design matrices via `Matrix::sparse.model.matrix()` and solve the
  least-squares problem with sparse QR or Cholesky methods. The
  reporter-by-perturbation design is highly structured and mostly zeros, so a
  dense `lm()` fit wastes memory.
- Use high-dimensional fixed-effect regression backends such as `fixest` or
  `lfe` to partial out reporter, perturbation, plate, replicate, concentration,
  or batch effects efficiently. These packages are designed for large factor
  structures. The main implementation question is how best to recover the
  reporter-by-perturbation contrasts and their covariance estimates in the
  notation used by DStressR.
- Use a Frisch-Waugh-Lovell strategy: residualize the response and the
  reporter-by-perturbation contrast columns with respect to technical and main
  effects, then solve only the reduced interaction problem. This would keep the
  mathematical estimand unchanged while avoiding a large dense joint design.
- For screens with many perturbations, keep the current reporter-wise scalable
  path as the default and reserve the full interaction model for cases where
  perturbation-level covariates, such as concentration, are scientifically
  central and the interaction block is still manageable.

A future implementation should therefore separate the statistical model from the
linear-algebra backend. The user-facing model could remain the same, while
`backend = "dense"`, `backend = "sparse"`, or `backend = "fixed_effects"`
selects the computational strategy. The default should continue to be the
scalable reporter-wise model for wide screens, with the interaction model
recommended for tall designs with a small number of perturbations or for
analyses that require explicit perturbation-level covariate adjustment.

For the PFIbox promoter-activity analysis, the intended dose-adjusted scalable
model is reporter-wise:

\[
r_i =
\eta_{a_i}
+
\theta_{a_i,\mathrm{drug}(i)}
+
\gamma_{a_i}\,\mathrm{dose}_i
+
\varepsilon_i .
\]

Thus each reporter/promoter has its own dose trend \(\gamma_a\). This is
scientifically more natural than a shared dose coefficient because fluorescence
or promoter activity can have reporter-specific concentration dependence. It
also means that the scalable reporter-wise implementation is the intended model,
not merely an approximation to a joint interaction model with a single shared
dose coefficient.

The practical issue is rank deficiency. If dose is completely determined by the
perturbation labels within a reporter, the perturbation indicators and dose
column are collinear. The model remains estimable only after dropping aliased
columns. The package should therefore use the fast closed-form reporter-wise
solver when the design is full rank and fall back to reporter-wise QR fitting
when the local design is rank deficient. This keeps the estimand reporter-wise
while avoiding the dense global interaction model.
