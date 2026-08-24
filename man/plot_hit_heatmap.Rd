\name{plot_hit_heatmap}
\alias{plot_hit_heatmap}
\title{Heatmap of significant DStressR hits}
\usage{
plot_hit_heatmap(
  table,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter",
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  show_compound_ids = TRUE,
  top_n_compounds = 160,
  fdr = 0.05,
  lfc = 0,
  drop_empty_compounds = TRUE,
  promoter_order = NULL,
  order_rows = c("global", "cluster", "input", "frequency"),
  order_cols = c("frequency", "cluster", "input"),
  clip_quantile = 0.98,
  color_limit = NULL,
  show_compound_labels = NULL,
  top_compound_labels = 40,
  compound_label_score = NULL,
  compound_label_min_gap = NULL,
  compound_label_angle = 45,
  title = "DStressR significant-hit matrix",
  subtitle = NULL,
  xlab = "Compounds",
  ylab = "Promoters",
  legend_title = effect,
  low = "#2166AC",
  mid = "white",
  high = "#B2182B"
)
}
\arguments{
\item{table}{A data frame with one row per promoter-compound pair.}

\item{effect}{Effect-size column shown by color for significant hits.}

\item{padj}{Adjusted p-value column used for hit calls.}

\item{promoter, compound}{Columns identifying promoters and compounds.}

\item{compound_label}{Optional human-readable compound-name column. Defaults to
\code{compound}.}

\item{show_compound_ids}{If \code{TRUE}, append compound IDs in square brackets
to compound labels.}

\item{top_n_compounds}{If finite, show only compounds with the strongest hit
evidence, ranked by hit count and effect size. Use \code{Inf} to show all
compounds.}

\item{fdr}{FDR threshold for hit highlighting.}

\item{lfc}{Minimum absolute effect size for hit highlighting.}

\item{drop_empty_compounds}{If \code{TRUE}, remove compounds with no
significant hits from the displayed matrix.}

\item{promoter_order}{Optional global promoter order used when
\code{order_rows = "global"}. If omitted, the option
\code{DStressR.promoter_order} is used when set; otherwise known DStressR paper
promoters are shown in their manuscript order and remaining promoters are
sorted alphabetically.}

\item{order_rows, order_cols}{Ordering strategy for promoters and compounds.
Use \code{"global"} for the package-wide promoter order, \code{"input"} to
preserve the input/factor order, \code{"frequency"} to order by number of hits,
or \code{"cluster"} for hierarchical clustering of the hit matrix.}

\item{clip_quantile}{Quantile of absolute significant effects used to clip the
color scale. Set to \code{1} to use the observed maximum.}

\item{color_limit}{Optional positive color-scale limit. If supplied, values are
clipped to \code{[-color_limit, color_limit]}; otherwise the limit is computed
from \code{clip_quantile}.}

\item{show_compound_labels}{If \code{TRUE}, draw all x-axis compound labels. If
\code{FALSE}, suppress x-axis compound labels. The default labels the
\code{top_compound_labels} compounds with largest absolute column sums, or all
compounds when fewer are plotted.}

\item{top_compound_labels}{Number of highest-signal compounds to label when
\code{show_compound_labels = NULL}.}

\item{compound_label_score}{Optional numeric score used to choose the
top-labelled compounds when \code{show_compound_labels = NULL}. If named,
values are matched to compound labels; otherwise the order must match the
displayed matrix columns.}

\item{compound_label_min_gap}{Minimum number of matrix columns between
automatically selected labels. The default chooses a gap from the displayed
matrix size and \code{top_compound_labels}.}

\item{compound_label_angle}{Angle used for visible compound labels.}

\item{title, subtitle, xlab, ylab}{Plot labels.}

\item{legend_title}{Colorbar title. Defaults to the selected \code{effect}
column.}

\item{low, mid, high}{Colors for negative, zero, and positive hit effects.}
}
\value{
A \code{ggplot} object with attributes \code{hit_matrix}, \code{hit_summary},
\code{plotted_pairs}, and \code{color_limit}.
}
\description{
Shows a promoter-by-compound effect matrix with significant pairs highlighted by
color and non-significant pairs shown as a light background. This plot is a
compact companion to \code{\link[=plot_response_heatmap]{plot_response_heatmap()}}
for inspecting the discovered hit structure.
}
