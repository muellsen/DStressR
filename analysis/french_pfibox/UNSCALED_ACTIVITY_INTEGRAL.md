# Unscaled activity integral

The retained manuscript analysis reconstructs unscaled activity from deposited
`constit.normalized` fluorescence for each colony. It does not use the
min-max-scaled `promoter.activity` values as its response.

Let q[k] denote a colony's biomass-normalized fluorescence, sampled every
dt = 1/12 hour. Retain differences q[k+1]-q[k] for k = 20,...,216,
then apply quadratic LOESS with span 0.5 on this 197-sample grid. Activity
per hour is the smoothed difference divided by dt. The response is its
signed trapezoidal integral over the retained interval, without min-max
scaling, logarithmic transformation, or clipping negative activity.

The default run integrates all 197 retained samples. Their right-endpoint
time coordinates are 100 minutes through 18 hours, assuming the public
five-minute grid starts at zero. This is the index alignment recovered
numerically, not an additional two-hour exclusion. Requested subwindows
can instead be selected with PFIBOX_TIME_WINDOW=requested and
PFIBOX_T_START_H / PFIBOX_T_END_H; endpoints use the available grid.

The integration is computed using the shared linear LOESS operator and
trapezoidal weights. The dt factors from differentiation and integration
cancel. Against the 940 individually reconstructed audit curves, maximum
absolute integral discrepancy was 2.19e-11 (floating-point arithmetic).

The fit uses the existing reporter-wise model with treatment levels and
the existing colony-replicate-position covariate, empirical-Bayes moderation,
and pooled untreated reference observations. Each source colony-condition
appears once; controls are not duplicated for each dose. The 13 untreated
source matrices are not described as independent biological replicates.
This run does not change the reference structure or add spatial adjustment.

All 1,797 named reporters and 33 non-reference treatment levels are retained.
There are 273,792 finite colony responses and 59,301 fitted total effects.
Constant upstream traces have a valid zero integral and are retained, even
when min-max scaling would have returned undefined deposited activity.
Thus coverage differs from the earlier scaled-activity runs.

Heatmap: unscaled total effects, no reporter-wise standardization; complete
linkage with correlation distance on both axes. The one constant reporter
profile is excluded from correlation clustering and retained in inference.
The current 1,796-by-33 displayed matrix is complete.
Dose-code labels are not converted into unverified concentrations.
The symmetric color range clips at the 98.5th percentile of absolute effects;
the saved diagnostic artifact retains the full unclipped matrix and row order.

This is a new DStressR response analysis, not a reproduction of the published
Figure 3 clustering. The French paper defines dynamic activity and its
min-max-scaled visualization; an exact published PA-integral analysis has
not been established. The manuscript now describes this integral response;
abandoned exploratory analyses have been removed.

## Figure 6 composite

`build_selected_activity_figure.R --target=total` saves the retained composite
under `selected_four_figure6/`. Panels show the total-effect matrix, four
colony-level time-course examples in a 2-by-2 layout, raw total-effect
p-values, and the corresponding volcano plot.
Significance uses global BH adjustment across all 59,301 tests (5,101 hits).
The volcano is capped at -log10(p) = 50; capped annotations give the true
-log10(p), as in the previous figure.

Correlation clustering excludes one constant effect profile, leaving 1,796
reporters; inference retains it. The heatmap color scale remains unstandardized.
The supplementary comparison uses major temporal-cluster labels from
Table S3. The generated `figure6_cluster_mapping_audit.tsv`
records the matched/unmatched counts and adjusted Rand index. Low ARI measures
partition disagreement, not statistical independence between biological traits.

```bash
Rscript analysis/french_pfibox/build_selected_activity_figure.R --target=total
```

## Selected time courses

`plot_integrated_activity_timecourses.R` shows individual unscaled PA curves
and pointwise condition medians. Untreated curves include the pooled reference
sources used by the fitted model; these medians are descriptive curves, not
covariate-adjusted predictions. Tests concern integrated responses.

For `PFIBOX_TIMECOURSE_SELECTION=directional_volcano`, selection is separately
performed among positive and negative effects with global BH q < 0.05. Absolute
effect size and -log10(raw p) are converted to within-direction percentile
ranks; their geometric mean defines a display score. The ten highest-scoring
eligible pairs are selected greedily, with at most one pair per reporter and
two pairs per antibiotic in each direction. This is an explicitly selected
illustration, not an additional inferential procedure. Pair tables retain the
component ranks, score, effects, p-values, and q-values.

```bash
PFIBOX_TIMECOURSE_SELECTION=directional_volcano \
Rscript analysis/french_pfibox/plot_integrated_activity_timecourses.R
```

```bash
PFIBOX_RESPONSE_SUMMARY=unscaled_integral \
PFIBOX_TIME_WINDOW=all_available \
PFIBOX_REPORTER_LIMIT=2000 PFIBOX_PERTURBATION_LIMIT=100 \
PFIBOX_DESIGN=replicate \
Rscript analysis/french_pfibox/run_promoter_activity_analysis.R

Rscript analysis/french_pfibox/prepare_integrated_activity_supplement.R
```
