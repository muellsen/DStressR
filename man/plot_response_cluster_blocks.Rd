\name{plot_response_cluster_blocks}
\alias{plot_response_cluster_blocks}
\title{Clustered block map of a DStressR promoter-by-compound response matrix}
\usage{
plot_response_cluster_blocks(
  table,
  value = "specific_effect",
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  show_compound_ids = TRUE,
  n_promoter_clusters = 6,
  n_compound_clusters = 14,
  missing_value = 0,
  clip_quantile = 0.98,
  show_counts = TRUE,
  title = "DStressR clustered response map",
  subtitle = NULL,
  xlab = "Compound clusters",
  ylab = "Promoter clusters",
  low = "#2166AC",
  mid = "white",
  high = "#B2182B"
)
}
\arguments{
\item{table}{A data frame with one row per promoter-compound pair.}

\item{value}{Numeric response/effect column to summarize.}

\item{promoter, compound}{Columns identifying promoters and compounds.}

\item{compound_label}{Optional human-readable compound-name column. Defaults to
\code{compound}.}

\item{show_compound_ids}{If \code{TRUE}, append compound IDs in square brackets
to compound labels before clustering.}

\item{n_promoter_clusters, n_compound_clusters}{Number of dendrogram clusters to
use for promoters and compounds.}

\item{missing_value}{Value used only for clustering missing matrix entries.
Block summaries are still computed from observed finite values.}

\item{clip_quantile}{Quantile of absolute block means used to clip the color
scale. Set to \code{1} to use the observed maximum.}

\item{show_counts}{If \code{TRUE}, annotate each tile with the number of
compounds in that compound cluster.}

\item{title, subtitle, xlab, ylab}{Plot labels.}

\item{low, mid, high}{Colors for negative, zero, and positive responses.}
}
\value{
A \code{ggplot} object with attributes \code{response_matrix},
\code{promoter_clusters}, \code{compound_clusters}, \code{block_summary},
\code{row_hclust}, and \code{col_hclust}.
}
\description{
Hierarchically clusters promoters and compounds, cuts the dendrograms into
interpretable groups, and plots the mean response for each promoter-cluster by
compound-cluster block. This is useful as a compact overview when the full
compound library is too large for individual compound labels.
}
