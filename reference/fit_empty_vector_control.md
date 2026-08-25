# Reproduce Empty Vector Control normalization

This function implements the Salmonella StressRegNet workflow in which
reporter activity is normalized against an Empty Vector Control (EVC)
reporter measured for the same perturbation/library well. It starts from
a long expression table, subtracts a perturbation-specific EVC average
from each reporter-replicate value, estimates reporter-replicate DMSO
null distributions, and applies the original conservative replicate
aggregation.

## Usage

``` r
fit_empty_vector_control(
  data,
  reporter = "reporter",
  perturbation = "srn_code",
  replicate = "replicate",
  response = "log2.lux.normed.centered",
  empty_vector_reporter = "PEVC3",
  control,
  exclude = character(),
  remove_reporters = character(),
  fdr = 0.05,
  require_complete_empty_vector = TRUE
)
```

## Arguments

- data:

  Long expression table with one row per reporter-perturbation-replicate
  observation.

- reporter, perturbation, replicate:

  Column names identifying reporter, perturbation/library well, and
  replicate.

- response:

  Column containing the expression value to normalize. For the
  Salmonella workflow this is `log2.lux.normed.centered`.

- empty_vector_reporter:

  Reporter/control strain used as the Empty Vector reference. The
  original Salmonella workflow uses `PEVC3`.

- control:

  Character vector of perturbation/library-well IDs used as DMSO
  controls for the null distribution.

- exclude:

  Character vector of perturbation/library-well IDs removed before
  normalization and hit calling, for example noisy DMSO wells.

- remove_reporters:

  Reporters removed before normalization, for example failed reporter
  strains.

- fdr:

  FDR threshold used to assign the `hit` class in the pair-level table.

- require_complete_empty_vector:

  If `TRUE`, require all EVC replicate values for a perturbation to be
  finite before computing the EVC average. This matches the original
  workflow's effective behavior with two PEVC3 replicates.

## Value

A list of class `destress_empty_vector` with `replicate_results`,
`pair_results`, `empty_vector_reference`, `control`, and `exclude`
components.
