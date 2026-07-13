\name{plot_response_clustered_heatmap}
\alias{plot_response_clustered_heatmap}
\title{Clustered heatmap with promoter and compound dendrograms}
\usage{
plot_response_clustered_heatmap(
  table,
  value = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  show_compound_ids = TRUE,
  top_n_compounds = 400,
  n_promoter_clusters = 6,
  n_compound_clusters = 14,
  missing_value = 0,
  clip_quantile = 0.98,
  file = NULL,
  width = 14,
  height = 8,
  res = 300,
  title = "DStressR clustered response heatmap",
  subtitle = NULL,
  show_rownames = TRUE,
  show_colnames = FALSE,
  low = "#2166AC",
  mid = "white",
  high = "#B2182B"
)
}
\arguments{
\item{table}{A data frame with one row per promoter-compound pair.}

\item{value}{Numeric response/effect column to show in the heatmap.}

\item{promoter, compound}{Columns identifying promoters and compounds.}

\item{compound_label}{Optional human-readable compound-name column. Defaults to
\code{compound}.}

\item{show_compound_ids}{If \code{TRUE}, append compound IDs in square brackets
to compound labels before clustering.}

\item{top_n_compounds}{If finite, show only the top compounds by mean absolute
response. Use \code{Inf} to show all compounds.}

\item{n_promoter_clusters, n_compound_clusters}{Number of dendrogram clusters
returned in the cluster assignment tables.}

\item{missing_value}{Value used only for clustering missing matrix entries.
Heatmap cells with missing values are left missing.}

\item{clip_quantile}{Quantile of absolute response values used to clip the
color scale. Set to \code{1} to use the observed maximum.}

\item{file}{Optional output file. Supports \code{.png} and \code{.pdf}. If
\code{NULL}, the plot is drawn on the active graphics device.}

\item{width, height}{Plot size in inches when \code{file} is supplied.}

\item{res}{PNG resolution in dots per inch.}

\item{title, subtitle}{Plot title and subtitle.}

\item{show_rownames, show_colnames}{Whether to draw row and column labels.}

\item{low, mid, high}{Colors for negative, zero, and positive responses.}
}
\value{
Invisibly returns a list containing the response matrix, clustering objects,
ordered matrix, cluster assignments, and color limit.
}
\description{
Draws a clustered promoter-by-compound response heatmap with hierarchical trees
on both axes. Unlike \code{\link{plot_response_cluster_blocks}}, this keeps the
individual matrix cells visible and uses the dendrograms to reveal structure
without collapsing the data into coarse blocks.
}
