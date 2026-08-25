# Summarize significant reporter-perturbation hits

Builds compact pair-, reporter-, and perturbation-level summaries from a
DStressR result table. Hits are defined by an adjusted p-value threshold
and, optionally, a minimum absolute effect size.

## Usage

``` r
summarize_hits(
  table,
  effect = "specific_effect",
  padj = "specific_padj_by_reporter",
  reporter = "reporter",
  perturbation = "perturbation",
  perturbation_label = perturbation,
  fdr = 0.05,
  lfc = 0
)
```

## Arguments

- table:

  A data frame with one row per reporter-perturbation pair.

- effect:

  Effect-size column used for hit direction and effect summaries.

- padj:

  Adjusted p-value column used for hit calls.

- reporter, perturbation:

  Columns identifying reporters and perturbations.

- perturbation_label:

  Optional human-readable perturbation-name column. Defaults to
  `perturbation`.

- fdr:

  FDR threshold for hit calls.

- lfc:

  Minimum absolute effect size for hit calls.

## Value

A list of class `destress_hit_summary` with pair-level hits and
reporter- and perturbation-level summaries.
