\name{fit_effect_mixture}
\alias{fit_effect_mixture}
\title{Fit a three-part empirical-null mixture to promoter-compound effects}
\usage{
fit_effect_mixture(
  table,
  value = "specific_effect",
  promoter = "promoter",
  df = 4,
  max_iter = 2000,
  tol = 1e-06,
  min_scale = 1e-04,
  min_prior = 1e-04,
  padj_method = "BH"
)
}
\arguments{
\item{table}{A data frame with one row per promoter-compound pair.}

\item{value}{Numeric effect column, usually \code{specific_effect} or a
centered DStressR EB effect column.}

\item{promoter}{Column identifying promoters.}

\item{df}{Degrees of freedom for each Student-t component. Smaller values give
heavier tails.}

\item{max_iter}{Maximum EM iterations per promoter.}

\item{tol}{Relative log-likelihood convergence tolerance.}

\item{min_scale}{Lower bound for component scale.}

\item{min_prior}{Lower bound for component mixing proportions.}

\item{padj_method}{Multiple-testing correction method passed to
\code{\link[stats:p.adjust]{stats::p.adjust()}}, applied within promoter to
empirical-null p-values.}
}
\value{
The input table with posterior probabilities, local FDR, posterior non-null
probability, local-FDR q-values within promoter, empirical-null p-values,
promoter-wise adjusted p-values, and posterior class appended. The
promoter-level fitted parameters are available as
\code{attr(result, "mixture_summary")}.
}
\description{
Fits, separately for each promoter, a three-component Student-t mixture to
adjusted promoter-compound effects. The ordered components are interpreted as
repressed, null, and activated effects. This second-stage model is intended for
empirical-null calibration after the first-stage DStressR model has already
adjusted growth, technical factors, compound-wide effects, and promoter-specific
variance.
}
