# Call differential stress hits

Call differential stress hits

## Usage

``` r
call_hits(
  table,
  fdr = 0.05,
  lfc = 0,
  effect = "specific_effect",
  padj = "specific_padj"
)
```

## Arguments

- table:

  Result table from
  [`results()`](https://bio-datascience.github.io/DStressR/reference/results.md).

- fdr:

  FDR threshold.

- lfc:

  Minimum absolute effect size.

- effect:

  Effect column, usually `specific_effect` or `total_effect`.

- padj:

  Adjusted p-value column.
