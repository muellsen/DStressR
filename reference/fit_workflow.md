# Fit a named DStressR workflow

`fit_workflow()` is a compatibility entry point for named DStressR
workflows. The `workflow` argument selects a
[`fit_destress()`](https://bio-datascience.github.io/DStressR/reference/fit_destress.md)
preset: `model` fits the model-based DStressR workflow, `median_polish`
reproduces the legacy median-polish workflow, and `empty_vector_control`
reproduces the empty-vector-control workflow.

## Usage

``` r
fit_workflow(
  data,
  workflow = c("model", "median_polish", "empty_vector_control"),
  ...
)
```

## Arguments

- data:

  For `workflow = "model"`, a `destress_assay` produced by
  [`prepare_assay()`](https://bio-datascience.github.io/DStressR/reference/prepare_assay.md).
  For the compatibility workflows, a long expression table.

- workflow:

  One of
  [`destress_workflows()`](https://bio-datascience.github.io/DStressR/reference/destress_workflows.md).
  Hyphenated names and common aliases such as `"destress"`,
  `"median-polish"`, and `"evc"` are accepted.

- ...:

  Arguments passed to the selected workflow engine.

## Value

The fitted workflow object returned by the selected engine, with a
`destress_workflow` attribute naming the workflow used.

## Details

New analyses should prefer
[`fit_destress()`](https://bio-datascience.github.io/DStressR/reference/fit_destress.md)
directly so the staged statistical choices can be made explicit.
