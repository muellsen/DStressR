#' List available DStressR workflows
#'
#' @return A character vector of workflow names accepted by [fit_workflow()].
#' @export
destress_workflows <- function() {
  c("model", "median_polish", "empty_vector_control")
}

normalize_workflow <- function(workflow) {
  if (missing(workflow) || length(workflow) != 1 || is.na(workflow) || !nzchar(workflow)) {
    stop("`workflow` must be one workflow name.", call. = FALSE)
  }
  workflow <- gsub("-", "_", tolower(workflow), fixed = TRUE)
  aliases <- c(
    destress = "model",
    model_based = "model",
    medianpolish = "median_polish",
    empty_vector = "empty_vector_control",
    evc = "empty_vector_control"
  )
  if (workflow %in% names(aliases)) {
    workflow <- aliases[[workflow]]
  }
  choices <- destress_workflows()
  if (!workflow %in% choices) {
    stop(
      "Unknown workflow `", workflow, "`. Available workflows are: ",
      paste(choices, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  workflow
}

#' Fit a named DStressR workflow
#'
#' `fit_workflow()` is a compatibility entry point for named DStressR
#' workflows. The `workflow` argument selects a [fit_destress()] preset:
#' `model` fits the model-based DStressR workflow, `median_polish` reproduces
#' the legacy median-polish workflow, and `empty_vector_control` reproduces the
#' empty-vector-control workflow.
#'
#' New analyses should prefer [fit_destress()] directly so the staged
#' statistical choices can be made explicit.
#'
#' @param data For `workflow = "model"`, a `destress_assay` produced by
#'   [prepare_assay()]. For the compatibility workflows, a long expression
#'   table.
#' @param workflow One of [destress_workflows()]. Hyphenated names and common
#'   aliases such as `"destress"`, `"median-polish"`, and `"evc"` are accepted.
#' @param ... Arguments passed to the selected workflow engine.
#' @return The fitted workflow object returned by the selected engine, with a
#'   `destress_workflow` attribute naming the workflow used.
#' @export
fit_workflow <- function(data,
                         workflow = c("model", "median_polish", "empty_vector_control"),
                         ...) {
  workflow <- normalize_workflow(workflow[1])
  fit <- fit_destress(data, preset = workflow, ...)
  attr(fit, "destress_workflow") <- workflow
  fit
}
