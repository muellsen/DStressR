\name{fit_median_polish}
\alias{fit_median_polish}
\title{Reproduce the legacy median-polish workflow}
\usage{
fit_median_polish(
  data,
  promoter = "promoter",
  compound = "srn_code",
  libplate = "libplate",
  replicate = "replicate",
  response = "log2.auc.16hmeasured.normed",
  control,
  exclude = character(),
  fdr = 0.05,
  maxiter = 1000,
  eps = 1e-08
)
}
\arguments{
\item{data}{Long expression table with one row per
promoter-compound-replicate observation.}

\item{promoter, compound, libplate, replicate}{Column names identifying the
promoter, compound/library well, library plate, and replicate.}

\item{response}{Column containing the already growth-normalized log2
response, for example \code{log2.auc.16hmeasured.normed}.}

\item{control}{Character vector of compound/library-well IDs used as DMSO
controls.}

\item{exclude}{Character vector of compound/library-well IDs to remove before
median polishing and hit calling, for example noisy DMSO wells.}

\item{fdr}{FDR threshold used to assign the \code{hit} class in the pair-level
table.}

\item{maxiter, eps}{Passed to \code{\link[stats:medpolish]{stats::medpolish()}}.}
}
\value{
A list of class \code{destress_median_polish} with \code{replicate_results},
\code{pair_results}, \code{polished_matrix}, and \code{medpolish} components.
}
\description{
This function implements the original median-polish hit-determination
workflow used for the Campylobacter promoter-library screen. It starts from a
long expression table, centers each promoter-library-plate-replicate group by
its DMSO wells, applies \code{\link[stats:medpolish]{stats::medpolish()}} to the resulting
promoter-libplate-replicate by compound matrix, and computes z-test p-values
from the polished DMSO residual distribution.
}
\details{
The promoter-compound hit table follows the original conservative replicate
aggregation: DMSO and excluded control wells are removed, the largest
replicate-level p-value is retained for each promoter-compound pair, and
p-values are BH-adjusted within promoter.
}
