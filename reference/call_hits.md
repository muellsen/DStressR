# Call differential stress hits

Call differential stress hits

## Usage

``` r
call_hits(
  table,
  fdr = 0.05,
  lfc = 0,
  effect = "total_effect",
  padj = sub("_effect$", "_padj", effect)
)
```

## Arguments

- table:

  Result table from
  [`results()`](https://muellsen.github.io/DStressR/reference/results.md).

- fdr:

  FDR threshold.

- lfc:

  Minimum absolute effect size.

- effect:

  Effect column. Defaults to `total_effect`; use `relative_effect` for
  deviations from the reporter-panel average.

- padj:

  Adjusted p-value column.
