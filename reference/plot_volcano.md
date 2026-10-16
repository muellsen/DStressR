# Volcano plot for DStressR reporter-perturbation hits

Creates a standard volcano plot from a DStressR result table. The x-axis
is a reporter-perturbation effect size and the y-axis is the negative
log10 adjusted p-value. Significant hits are emphasized and the most
significant reporter-perturbation pairs can be annotated.

## Usage

``` r
plot_volcano(
  table,
  effect = "total_effect",
  padj = sub("_effect$", "_padj", effect),
  pvalue = NULL,
  y = c("auto", "pvalue", "padj"),
  reporter = "reporter",
  perturbation = "perturbation",
  perturbation_label = perturbation,
  fdr = 0.05,
  lfc = 0,
  top_n = 12,
  top_reporters = 6,
  color_by = c("significance", "reporter"),
  title = "DStressR volcano plot",
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  label_by = c("pair", "reporter", "perturbation"),
  label_hits_only = TRUE,
  max_label_chars = 46,
  repel_labels = TRUE,
  point_size = 1.35,
  point_alpha = 0.78,
  hit_color = "#D95F02",
  background_color = "#D1D5DB",
  call_labels = c("Not called", "Called")
)
```

## Arguments

- table:

  A data frame with one row per reporter-perturbation pair.

- effect:

  Effect-size column to plot on the x-axis.

- padj:

  Adjusted p-value column used for FDR calls.

- pvalue:

  Optional raw p-value column. When supplied, the default y-axis is
  `-log10(pvalue)` while calls are still defined by `padj`.

- y:

  Which p-value scale to show on the y-axis. The default `"auto"` uses
  `pvalue` when available and otherwise uses `padj`.

- reporter, perturbation:

  Columns identifying the reporter and perturbation.

- perturbation_label:

  Optional column with human-readable perturbation names used for
  annotations. Defaults to `perturbation`.

- fdr:

  FDR threshold for hit highlighting.

- lfc:

  Minimum absolute effect size for hit highlighting.

- top_n:

  Number of pairs to annotate. By default these are selected from
  significant hits.

- top_reporters:

  Number of reporter groups to color. Remaining reporters are shown in
  grey when `color_by = "reporter"`.

- color_by:

  Color points by `"significance"` or by `"reporter"`.

- title, subtitle:

  Plot title and subtitle.

- xlab, ylab:

  Axis labels. Defaults to readable labels based on the selected
  columns.

- label_by:

  Label style for annotated points. The default, `"pair"`, labels top
  hits as reporter-perturbation pairs.

- label_hits_only:

  If `TRUE`, annotations are selected only from pairs passing the FDR
  and effect-size thresholds.

- max_label_chars:

  Maximum characters per annotation label. Longer labels are truncated
  with `...`. Use `Inf` to keep full labels.

- repel_labels:

  If `TRUE` and the optional `ggrepel` package is installed, use
  repelled labels for the annotated top hits.

- point_size, point_alpha:

  Point size and transparency.

- hit_color, background_color:

  Colors used when `color_by = "significance"`.

- call_labels:

  Two labels for non-called and called pairs when
  `color_by = "significance"`.

## Value

A `ggplot` object.

## Details

The defaults work with
[`results()`](https://muellsen.github.io/DStressR/reference/results.md)
followed by
[`adjust_pvalues()`](https://muellsen.github.io/DStressR/reference/adjust_pvalues.md)
or
[`call_hits()`](https://muellsen.github.io/DStressR/reference/call_hits.md).
For workflow comparison tables, pass the corresponding column names, for
example `effect = "destress_eb_effect_centered"`,
`padj = "estimated_alpha_eb_padj_by_reporter"`,
`perturbation = "srn_code"`, and `perturbation_label = "ProductName"`.
