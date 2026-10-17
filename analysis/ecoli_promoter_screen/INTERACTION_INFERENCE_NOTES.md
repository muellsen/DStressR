# Interaction-model inference

## Total Effects Used in Supplementary Note F

The comparison uses compound-versus-water contrasts of the full fitted mean,
including the compound main coefficient and its reporter interaction term.
Standard errors come from the full coefficient covariance. Total-effect tests
now honor `empirical_bayes = TRUE`, with moderation applied once to the raw
contrast variances and the common full-model residual degrees of freedom.
Rank-adjusted total tests reuse the moderated total uncertainty.

The saved result artifact supplies the raw contrasts needed for the existing
paper analysis. The preparer validates the earlier ordinary p-values before
reconstructing moderated total tests. Future analysis runs also save fit objects.

The pre-change reporter-relative-effect outputs and manuscript are preserved locally in
`analysis/outputs/inference_audit_2026-10-16/pre_interaction_total_effects.tar.gz`.

## Reporter-relative-Effect Covariance

The full-interaction path centers each treatment-versus-reference design
contrast across the selected reporters, then evaluates its variance using the
full coefficient covariance. The raw centered variance is stored as
`relative_var` and moderated once with the common model residual degrees of
freedom. This replaces the older reuse of uncentered total-effect uncertainty.
Tests compare the result with explicitly constructed joint contrasts in an
unbalanced design and cover moderation, rank adjustment, and a single-reporter
selection. Supplementary Note F uses total effects only and is unchanged by
this reporter-relative-effect correction. Earlier archived reporter-relative-effect outputs must
be regenerated before reuse.
