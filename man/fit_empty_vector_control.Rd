\name{fit_empty_vector_control}
\alias{fit_empty_vector_control}
\title{Reproduce Empty Vector Control normalization}
\usage{
fit_empty_vector_control(
  data,
  promoter = "promoter",
  compound = "srn_code",
  replicate = "replicate",
  response = "log2.lux.normed.centered",
  empty_vector_promoter = "PEVC3",
  control,
  exclude = character(),
  remove_promoters = character(),
  fdr = 0.05,
  require_complete_empty_vector = TRUE
)
}
\arguments{
\item{data}{Long expression table with one row per
promoter-compound-replicate observation.}

\item{promoter, compound, replicate}{Column names identifying promoter,
compound/library well, and replicate.}

\item{response}{Column containing the expression value to normalize. For the
Salmonella workflow this is \code{log2.lux.normed.centered}.}

\item{empty_vector_promoter}{Promoter/control strain used as the Empty Vector
reference. The original Salmonella workflow uses \code{PEVC3}.}

\item{control}{Character vector of compound/library-well IDs used as DMSO
controls for the null distribution.}

\item{exclude}{Character vector of compound/library-well IDs removed before
normalization and hit calling, for example noisy DMSO wells.}

\item{remove_promoters}{Promoters removed before normalization, for example
failed reporter strains.}

\item{fdr}{FDR threshold used to assign the \code{hit} class in the pair-level
table.}

\item{require_complete_empty_vector}{If \code{TRUE}, require all EVC replicate
values for a compound to be finite before computing the EVC average. This
matches the original workflow's effective behavior with two PEVC3
replicates.}
}
\value{
A list of class \code{destress_empty_vector} with \code{replicate_results},
\code{pair_results}, \code{empty_vector_reference}, \code{control}, and
\code{exclude} components.
}
\description{
This function implements the Salmonella StressRegNet workflow in which
promoter activity is normalized against an Empty Vector Control (EVC) reporter
measured for the same compound/library well. It starts from a long expression
table, subtracts a compound-specific EVC average from each
promoter-replicate value, estimates promoter-replicate DMSO null
distributions, and applies the original conservative replicate aggregation.
}
