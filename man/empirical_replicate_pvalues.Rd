\name{empirical_replicate_pvalues}
\alias{empirical_replicate_pvalues}
\title{Empirical p-values from replicate-averaged compound effects and DMSO nulls}
\usage{
empirical_replicate_pvalues(
  table,
  value,
  promoter = "promoter",
  compound = "compound",
  control,
  replicate = NULL,
  strata = NULL,
  min_replicates = 2,
  min_null = 5,
  permutation = FALSE,
  B = 1000,
  seed = NULL,
  alternative = c("two.sided", "greater", "less"),
  padj_method = "BH"
)
}
\arguments{
\item{table}{Replicate-level data frame.}

\item{value}{Numeric value column to test, for example a DStressR adjusted
effect such as \code{destress_eb_effect_centered}.}

\item{promoter, compound}{Columns identifying promoters and compounds.}

\item{control}{Character vector of control compound IDs, usually DMSO wells.}

\item{replicate}{Optional replicate column. The function does not require
replicate labels for averaging, but the argument documents the intended
replicate-level input and is checked when supplied.}

\item{strata}{Optional columns defining matched null strata. Use
\code{strata = "libplate"} to compare compounds only to DMSO wells from the
same library plate.}

\item{min_replicates}{Minimum finite replicate values required for an averaged
compound or DMSO control value.}

\item{min_null}{Minimum number of averaged DMSO controls required to compute an
empirical p-value.}

\item{permutation}{If \code{TRUE}, compute an additional permutation p-value by
sampling replicate-sized DMSO sets within each matched null stratum.}

\item{B}{Number of permutation draws when \code{permutation = TRUE}.}

\item{seed}{Optional random seed for reproducible permutation p-values.}

\item{alternative}{One of \code{"two.sided"}, \code{"greater"}, or
\code{"less"}.}

\item{padj_method}{Multiple-testing correction method passed to
\code{\link[stats:p.adjust]{stats::p.adjust()}}, applied within promoter.}
}
\value{
A data frame with one row per non-control promoter-compound-stratum average,
including empirical p-values and promoter-wise adjusted p-values.
}
\description{
Computes empirical p-values by first averaging replicate-level values for each
promoter-compound-stratum combination, then comparing each averaged non-control
compound value to the corresponding distribution of averaged DMSO control values
from the same promoter and stratum, such as library plate.
Optionally, it also computes Monte Carlo permutation p-values by repeatedly
drawing replicate-sized DMSO sets from the same matched null stratum.
}
