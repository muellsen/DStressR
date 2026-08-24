\name{summarize_hits}
\alias{summarize_hits}
\title{Summarize significant promoter-compound hits}
\usage{
summarize_hits(
  table,
  effect = "specific_effect",
  padj = "specific_padj_by_promoter",
  promoter = "promoter",
  compound = "compound",
  compound_label = compound,
  fdr = 0.05,
  lfc = 0
)
}
\arguments{
\item{table}{A data frame with one row per promoter-compound pair.}

\item{effect}{Effect-size column used for hit direction and effect summaries.}

\item{padj}{Adjusted p-value column used for hit calls.}

\item{promoter, compound}{Columns identifying promoters and compounds.}

\item{compound_label}{Optional human-readable compound-name column. Defaults to
\code{compound}.}

\item{fdr}{FDR threshold for hit calls.}

\item{lfc}{Minimum absolute effect size for hit calls.}
}
\value{
A list of class \code{destress_hit_summary} with pair-level hits and
promoter- and compound-level summaries.
}
\description{
Builds compact pair-, promoter-, and compound-level summaries from a DStressR
result table. Hits are defined by an adjusted p-value threshold and, optionally,
a minimum absolute effect size.
}
