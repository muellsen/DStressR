# Verify public table identities and optionally refit without replacing paper results.
source(file.path("analysis", "_helpers.R"))
load_destress_package()
out <- analysis_output_dir("french_pfibox", "temporal_unscaled_integral_analysis")
public <- analysis_path("data-raw", "french_pfibox", "derived")
tag <- "reporters2000_perturbations100_Tall_available_treatment_replicate_scalable_temporal_unscaled_integral"
design <- read.delim(file.path(public, "experimental_design.tsv.gz"))
scalar <- read.delim(file.path(public, "scalar_responses.tsv.gz"))
stopifnot(nrow(design) == 273792L, !anyDuplicated(design$observation_id),
  !anyDuplicated(scalar$observation_id), setequal(design$observation_id, scalar$observation_id))
input <- cbind(design, scalar[match(design$observation_id, scalar$observation_id),
  setdiff(names(scalar), "observation_id")])
original <- read.delim(file.path(out, paste0("pfibox_pa_response_table_", tag, ".tsv")))
key <- paste(original$treatment, original$colony_position, sep = ":")
stopifnot(identical(input$observation_id, key),
  isTRUE(all.equal(input$value, original$y_pa, tolerance = 1e-12)),
  identical(input$Promoter, original$reporter))
message("Public metadata and scalar tables reproduce every saved colony response.")
if ("--refit" %in% commandArgs(trailingOnly = TRUE)) {
  assay <- prepare_assay(input, reporter = "Promoter", perturbation = "Perturbation",
    control = "Control", response = "value", replicate = "repeat_position")
  fit <- fit_destress(assay, technical = "repeat_position", testing = "moderated_t",
    adjustment = "global", interaction = FALSE, verbose = TRUE)
  fresh <- results(fit)
  saved <- read.delim(file.path(out, paste0("pfibox_pa_destress_pair_results_", tag, ".tsv")))
  idx <- match(paste(saved$reporter, saved$perturbation), paste(fresh$reporter, fresh$perturbation))
  stopifnot(!anyNA(idx))
  fresh <- fresh[idx, ]
  stopifnot(isTRUE(all.equal(saved$total_effect, fresh$total_effect, tolerance = 1e-10)),
    isTRUE(all.equal(saved$total_pvalue, fresh$total_pvalue, tolerance = 1e-10)),
    identical(saved$total_padj_global < .05, fresh$total_padj_global < .05))
  message("Public-table refit reproduces total effects, p-values and all 5,101 global BH calls.")
}
