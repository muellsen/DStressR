# Volcano plot for DStressR promoter-compound hits

Creates a standard volcano plot from a DStressR result table. The x-axis
is a promoter-compound effect size and the y-axis is the negative log10
adjusted p-value. Significant hits are emphasized, top promoter groups
can be colored, and the most significant promoter-compound pairs are
annotated.

## Usage

``` r
plot_volcano(
  table,
  effect = "specific_effect",
  padj = "specific_padj",
  pvalue = NULL,
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  fdr = 0.05,
  lfc = 0,
  top_n = 12,
  top_promoters = 6,
  title = "DStressR volcano plot",
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  label_by = c("pair", "promoter", "compound"),
  max_label_chars = 46,
  repel_labels = TRUE,
  point_alpha = 0.65
)
```

## Arguments

- table:

  A data frame with one row per promoter-compound pair.

- effect:

  Effect-size column to plot on the x-axis.

- padj:

  Adjusted p-value column to plot on the y-axis.

- pvalue:

  Optional raw p-value column used only if `padj = NULL`.

- promoter, compound:

  Columns identifying the promoter and compound.

- compound_label:

  Optional column with human-readable compound names used for
  annotations. Defaults to `compound`.

- fdr:

  FDR threshold for hit highlighting.

- lfc:

  Minimum absolute effect size for hit highlighting.

- top_n:

  Number of significant pairs to annotate.

- top_promoters:

  Number of promoter groups to color. Remaining promoters are shown in
  grey.

- title, subtitle:

  Plot title and subtitle.

- xlab, ylab:

  Axis labels. Defaults to readable labels based on the selected
  columns.

- label_by:

  Label style for annotated points. The default, `"pair"`, labels top
  hits as promoter-compound pairs.

- max_label_chars:

  Maximum characters per annotation label. Longer labels are truncated
  with `...`. Use `Inf` to keep full labels.

- repel_labels:

  If `TRUE` and the optional `ggrepel` package is installed, use
  repelled labels for the annotated top hits.

- point_alpha:

  Point transparency.

## Value

A `ggplot` object.

## Details

The defaults work with
[`results()`](https://bio-datascience.github.io/DStressR/reference/results.md)
followed by
[`adjust_pvalues()`](https://bio-datascience.github.io/DStressR/reference/adjust_pvalues.md)
or
[`call_hits()`](https://bio-datascience.github.io/DStressR/reference/call_hits.md).
For workflow comparison tables, pass the corresponding column names, for
example `effect = "destress_eb_effect_centered"`,
`padj = "estimated_alpha_eb_padj_by_promoter"`, `compound = "srn_code"`,
and `compound_label = "ProductName"`.
