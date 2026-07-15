# Adjust p-values within promoter

Adjust p-values within promoter

## Usage

``` r
adjust_pvalues(
  table,
  pvalue = "specific_pvalue",
  output = "specific_padj_by_promoter",
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
