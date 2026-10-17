# PFIbox source-stage and activity reconstruction audit

The local `gfpData.Rdata` is the Mendeley Data v1 object:
<https://data.mendeley.com/datasets/fmv6vv4bsf/1>.
Size: 1,311,202,733 bytes. SHA-256:
`41b3f09bc7a8c3c435af9d03d0e33cede018ba7a5fd0ec900cf6f7d07079ed63`.
The checksum and size match the official repository metadata.

The deposited `promoter.activity` trajectories are min-max-scaled per colony.
Amplitude for the integrated response is recovered from `constit.normalized`,
the deposited background-corrected, biomass-normalized fluorescence.

## Exact numerical recovery

`audit_french_figure3_inputs.R` fixes the two highest-membership promoters in
each of the twelve published clusters before inspecting trajectories. These
24 reporters yield 940 complete activity curves across the thirteen untreated
matrices. Checked colony-position identities match the public plate layout.

`audit_activity_reconstruction.R` compares declared shared processing rules:

1. Difference the 217-sample normalized fluorescence trace.
2. Retain difference indices 20--216 (197 samples).
3. Apply quadratic LOESS with span 0.5 on this retained grid.
4. Min-max-scale the fitted values for comparison to deposited activity.

This reproduces every checked activity trajectory to floating-point precision:
median RMSE 2.61e-16, maximum RMSE 1.41e-15. Smoothing before trimming does
not reproduce the stored curves. The primary DStressR response integrates the
unscaled activity; scaling is used only for this recovery check.

The smoothed differences become activity per hour after division by the
five-minute interval, dt = 1/12 hour. Trapezoidal integration multiplies by dt,
so these two factors cancel. The right-endpoint grid is 100 minutes through
18 hours under the source's five-minute acquisition convention.

The retained audit cache, all tested rules, and the selected-colony identifiers
are in `analysis/outputs/french_pfibox/figure3_input_audit/`.

## Published cluster labels

French et al.'s Table S3 supplies 1,481 promoter major-cluster assignments and
twelve fuzzy membership values for untreated temporal profiles. Our mapping
uses these reported assignments directly. It does not claim to reproduce the
historical Figure 3 colony selection, aggregation or Mfuzz input trajectories.
The current effect-profile comparison has 1,478 matched promoters; its low
adjusted Rand index is a partition comparison, not a test of independence.

## Processing-code provenance

The public scripts at <https://github.com/sfrench007/pfibox> explain min-max
scaling and derivative smoothing. They are a later protocol example, including
MOPS sample data, not a frozen 2018 Figure 3 implementation. The exact recovery
above is checked numerically against the deposited multi-treatment object.
Abandoned attempts to reproduce the published figures have been removed.
