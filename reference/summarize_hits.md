# Summarize significant promoter-compound hits

Builds compact pair-, promoter-, and compound-level summaries from a
DStressR result table. Hits are defined by an adjusted p-value threshold
and, optionally, a minimum absolute effect size.

## Usage

``` r
summarize_hits(
  table,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter",
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  fdr = 0.05,
  lfc = 0
)
```

## Arguments

- table:

  A data frame with one row per promoter-compound pair.

- effect:

  Effect-size column used for hit direction and effect summaries.

- padj:

  Adjusted p-value column used for hit calls.

- promoter, compound:

  Columns identifying promoters and compounds.

- compound_label:

  Optional human-readable compound-name column. Defaults to `compound`.

- fdr:

  FDR threshold for hit calls.

- lfc:

  Minimum absolute effect size for hit calls.

## Value

A list of class `destress_hit_summary` with pair-level hits and
promoter- and compound-level summaries.
