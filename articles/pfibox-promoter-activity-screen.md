# PFIbox promoter-activity screen

This article documents the DStressR reanalysis of the high-density
PFIbox promoter-reporter screen of French, Coutts, and Brown (2018),
*Cell Systems*, <doi:10.1016/j.cels.2018.07.004>. The deposited data
contain colony-level time-course measurements for an *E. coli*
promoter-GFP library under antibiotic treatments. The manuscript
analysis reconstructs unscaled promoter activity from the deposited
biomass-normalized fluorescence, integrates each colony trajectory, and
applies DStressR effect estimation and inference.

## Data provenance

The source is [Mendeley Data version
1](https://data.mendeley.com/datasets/fmv6vv4bsf/1), DOI
10.17632/fmv6vv4bsf.1, under CC BY 4.0. The full object is approximately
1.31 GB. Download and checksum instructions are in the data-folder
README. The public processing-code repository is:

``` text
https://github.com/sfrench007/pfibox
```

The local workflow expects the deposited multi-treatment `gfpData.Rdata`
object under:

``` text
data-raw/french_pfibox/source/gfpData.Rdata
```

The French et al. time-course cluster assignment table can be placed in
the same source folder when reproducing the cluster-mapping figure.
Large source files and derived tables are intentionally not committed.
The provenance note is kept in:

``` text
data-raw/french_pfibox/README.md
```

## Response construction

The deposited `constit.normalized` matrices contain biomass-normalized
fluorescence q(t). The deposited `promoter.activity` matrices have
already been min-max-scaled per colony. We recover unscaled activity by
differencing q(t), retaining difference indices 20–216, and applying
quadratic LOESS with span 0.5 on the retained 197-sample grid. Scaling
these reconstructed curves recovers the deposited activity to
floating-point precision in the audit. The DStressR input is the signed
trapezoidal integral of the unscaled activity:

``` text
y_i = integral from t_start to T of reconstructed_PA_i(t) dt
```

The manuscript uses the complete retained grid: 100 minutes through 18
hours under the five-minute acquisition convention. No additional
two-hour cutoff is applied. The integration endpoint remains
configurable. The DStressR model uses treatment-level perturbations,
where each antibiotic-dose combination is a separate perturbation.
Repeated colony positions are included as reporter-wise adjustment terms
when supported by the filtered data.

## Reproduce the analysis

First audit the downloaded source object:

``` sh
Rscript analysis/french_pfibox/audit_pfibox_data.R
```

Then run the promoter-activity subset analysis. The values below
reproduce the large subset used for the manuscript figure.

``` sh
PFIBOX_REPORTER_LIMIT=2000 \
PFIBOX_PERTURBATION_LIMIT=100 \
PFIBOX_RESPONSE_SUMMARY=unscaled_integral \
PFIBOX_TIME_WINDOW=all_available \
PFIBOX_PERTURBATION_TARGET=treatment \
PFIBOX_DESIGN=replicate \
PFIBOX_MODEL_PATH=scalable \
Rscript analysis/french_pfibox/run_promoter_activity_analysis.R
```

This writes the supplied response table, fitted DStressR pair results,
and diagnostic summaries under:

``` text
analysis/outputs/french_pfibox/temporal_unscaled_integral_analysis/
```

Prepare ordinary tests and the 10-positive plus 10-negative illustrative
selection, then compute the clustering with the settings below:

``` sh
Rscript analysis/french_pfibox/prepare_integrated_activity_supplement.R
PFIBOX_TIMECOURSE_SELECTION=directional_volcano \
Rscript analysis/french_pfibox/plot_integrated_activity_timecourses.R

export PFIBOX_ANALYSIS_DIR=analysis/outputs/french_pfibox/temporal_unscaled_integral_analysis
export PFIBOX_RUN_TAG=reporters2000_perturbations100_Tall_available_treatment_replicate_scalable_temporal_unscaled_integral
export PFIBOX_CLUSTER_TOP_REPORTERS=2000
export PFIBOX_CLUSTER_LINKAGE=complete
export PFIBOX_CLUSTER_DISTANCE=correlation
export PFIBOX_CLUSTER_OUTPUT_SUFFIX=total_effect_correlation_complete
export PFIBOX_ADJUSTMENT_SCOPE=global
Rscript analysis/french_pfibox/compute_paper_diagnostics.R
Rscript analysis/french_pfibox/map_french_timecourse_clusters_to_destress.R
Rscript analysis/french_pfibox/build_selected_activity_figure.R --target=total
```

The moderation diagnostic in the supplement compares the same fitted
model under moderated and ordinary standard errors:

``` sh
Rscript analysis/french_pfibox/plot_pfibox_moderation_diagnostic.R --target=total
```

## DStressR model call

The derived tables under `data-raw/french_pfibox/derived/` provide
observation metadata and long-format scalar responses without loading
the large source. Their README records units, transformations,
attribution and checksums. They are excluded from the CRAN package
build. After joining by observation ID:

``` r

design <- read.delim("data-raw/french_pfibox/derived/experimental_design.tsv.gz")
scalar <- read.delim("data-raw/french_pfibox/derived/scalar_responses.tsv.gz")
input <- merge(design, scalar, by = "observation_id", sort = FALSE)

assay <- prepare_assay(
  input,
  reporter = "Promoter",
  perturbation = "Perturbation",
  control = "Control",
  response = "value",
  replicate = "repeat_position"
)

fit <- fit_destress(
  assay,
  technical = "repeat_position",
  testing = "moderated_t",
  adjustment = "global",
  interaction = FALSE,
  verbose = TRUE
)

tab <- results(fit, verbose = TRUE)
hits <- call_hits(
  tab,
  fdr = 0.05,
  effect = "total_effect",
  padj = "total_padj_global"
)
```

The fitted object stores the reporter-wise total-effect fits. The
[`results()`](https://muellsen.github.io/DStressR/reference/results.md)
step extracts total effects, reporter-relative effects, moderated
statistics, and adjusted p-values.

Both BH families are available in the results; the manuscript uses the
global family across all 59,301 tests, yielding 5,101 calls at 0.05. It
contains 273,792 finite colony responses from 1,797 reporters and 33
treatment levels. The thirteen untreated source arrays are pooled as
references, with each source colony-condition represented once.

Figure 6 uses correlation distance and complete linkage for the 1,796
nonconstant reporter profiles. Four manually chosen examples from the
saved 20-pair display selection are marked consistently in the heatmap,
shown in the 2-by-2 trajectory panel, and labelled in the volcano. They
are ompC / Ampicillin-18, msrB / Sulfamethoxazole-18, recN /
Ciprofloxacin-18 and rmf / NalidixicAcid-18. These illustrative choices
do not change the hypothesis tests.
