\name{plot_effect_histogram}
\alias{plot_effect_histogram}
\title{Histogram of DStressR promoter-compound effects}
\usage{
plot_effect_histogram(
  table,
  value = "specific_effect",
  promoter = "promoter",
  by = c("all", "promoter"),
  bins = 80,
  xlim = NULL,
  scales = "fixed",
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = "Count",
  fill = "#4E79A7",
  border = "white"
)
}
\arguments{
\item{table}{A data frame with one row per promoter-compound pair.}

\item{value}{Numeric effect column to plot.}

\item{promoter}{Column identifying promoters, used when \code{by = "promoter"}.}

\item{by}{Plot one pooled histogram (\code{"all"}) or promoter-faceted
histograms (\code{"promoter"}).}

\item{bins}{Number of histogram bins.}

\item{xlim}{Optional two-element x-axis limit.}

\item{scales}{Facet scale behavior for \code{by = "promoter"}.}

\item{title, subtitle, xlab, ylab}{Plot labels.}

\item{fill, border}{Histogram fill and border colors.}
}
\value{
A \code{ggplot} object.
}
\description{
Shows the empirical distribution of normalized promoter-compound effects, either
over all matrix entries or faceted by promoter.
}
