# Reproduce Empty Vector Control normalization

This function implements the Salmonella StressRegNet workflow in which
promoter activity is normalized against an Empty Vector Control (EVC)
reporter measured for the same compound/library well. It starts from a
long expression table, subtracts a compound-specific EVC average from
each promoter-replicate value, estimates promoter-replicate DMSO null
distributions, and applies the original conservative replicate
aggregation.

## Usage

``` r
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
```

## Arguments

- data:

  Long expression table with one row per promoter-compound-replicate
  observation.

- promoter, compound, replicate:

  Column names identifying promoter, compound/library well, and
  replicate.

- response:

  Column containing the expression value to normalize. For the
  Salmonella workflow this is `log2.lux.normed.centered`.

- empty_vector_promoter:

  Promoter/control strain used as the Empty Vector reference. The
  original Salmonella workflow uses `PEVC3`.

- control:

  Character vector of compound/library-well IDs used as DMSO controls
  for the null distribution.

- exclude:

  Character vector of compound/library-well IDs removed before
  normalization and hit calling, for example noisy DMSO wells.

- remove_promoters:

  Promoters removed before normalization, for example failed reporter
  strains.

- fdr:

  FDR threshold used to assign the `hit` class in the pair-level table.

- require_complete_empty_vector:

  If `TRUE`, require all EVC replicate values for a compound to be
  finite before computing the EVC average. This matches the original
  workflow's effective behavior with two PEVC3 replicates.

## Value

A list of class `destress_empty_vector` with `replicate_results`,
`pair_results`, `empty_vector_reference`, `control`, and `exclude`
components.
