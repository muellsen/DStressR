# Adjust p-values within reporter

Adjust p-values within reporter

## Usage

``` r
adjust_pvalues(
  table,
  pvalue = "total_pvalue",
  output = sub("_pvalue$", "_padj_by_reporter", pvalue),
  method = "BH"
)
```

## Arguments

- table:

  Result table from
  [`results()`](https://muellsen.github.io/DStressR/reference/results.md).

- pvalue:

  P-value column.

- output:

  Name of adjusted p-value column.

- method:

  Passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).
