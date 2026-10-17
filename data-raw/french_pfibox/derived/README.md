# PFIbox derived colony-level inputs

These compressed tables are derived from Shawn French's *PFIbox Data Set*,
Mendeley Data version 1, DOI <https://doi.org/10.17632/fmv6vv4bsf.1>,
released under **CC BY 4.0**. Cite the dataset and French, Coutts and Brown
(2018), *Cell Systems*, DOI <https://doi.org/10.1016/j.cels.2018.07.004>.
These files are transformed DStressR analysis inputs, not unchanged source data.
The source licence applies with this attribution and transformation notice.

## Files and schema

- `experimental_design.tsv.gz`: one row per retained source colony-condition.
  `Promoter` and the six `Plate`/`Row`/`Column` layout fields retain source
  names. `source_treatment` retains the original condition-entry name.
  `Perturbation` is the analysis coding: all thirteen `*-0` conditions map to
  `Control`, and each non-reference antibiotic-dose entry remains distinct.
  `source_dose_code` is an original label, not a claimed physical concentration.
  `repeat_position` is the within-promoter colony-position index used by the
  model, not an assertion of independent biological replication.
- `scalar_responses.tsv.gz`: long-format responses linked by `observation_id`.
  `summary` is `signed_pa_integral`, `value` is the supplied DStressR response,
  and the window and timepoint count are explicit. Activity is a derivative per
  hour of biomass-normalized fluorescence; integrating it returns the units
  of that normalized fluorescence, not calibrated mRNA or protein abundance.
- `manifest.tsv`: row counts, compressed sizes, derived-file MD5 checksums,
  source DOI, source SHA-256 and licence.

Both tables have 273,792 rows. The identifier combines the original treatment
name and original colony-position index. Each source observation occurs once.
Joining the two tables does not duplicate untreated controls across doses.

## Transformations

Starting from deposited `constit.normalized` fluorescence, take consecutive
differences, retain difference indices 20--216, and apply quadratic LOESS with
span 0.5 on that 197-sample grid. Convert the smoothed difference to a derivative
per hour using five-minute spacing. The response is its signed trapezoidal
integral over right-endpoint times 100 minutes through 18 hours. There is no
per-colony min-max scaling, logarithmic transform, or clipping of signed activity.
Rows require all retained source differences to be finite. Empty and designated
control-reporter labels are excluded as described in the analysis README.

## Read and model

```r
design <- read.delim("data-raw/french_pfibox/derived/experimental_design.tsv.gz")
scalar <- read.delim("data-raw/french_pfibox/derived/scalar_responses.tsv.gz")
input <- merge(design, scalar, by = "observation_id", sort = FALSE)
stopifnot(nrow(input) == 273792L)

assay <- DStressR::prepare_assay(
  input, reporter = "Promoter", perturbation = "Perturbation",
  control = "Control", response = "value", replicate = "repeat_position"
)
fit <- DStressR::fit_destress(
  assay, technical = "repeat_position", testing = "moderated_t",
  adjustment = "global", interaction = FALSE
)
tab <- DStressR::results(fit)
hits <- DStressR::call_hits(tab, effect = "total_effect",
  padj = "total_padj_global", fdr = 0.05)
```

Regenerate with `analysis/french_pfibox/export_public_response_tables.R` after
the retained response-construction workflow. The 1.31-GB original stays in the
Git-ignored source directory; these compact tables are intended for GitHub
and are excluded from the CRAN package build by the `data-raw/` build rule.
