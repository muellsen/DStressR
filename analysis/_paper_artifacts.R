# Lightweight contract for local paper outputs; this is not a package API.
write_paper_artifact <- function(path, results = list(), assays = list(),
                                 parameters = list(), diagnostics = list(),
                                 annotations = list(), provenance = list()) {
  for (table in results) {
    if (!is.data.frame(table) ||
        !all(c("reporter", "perturbation") %in% names(table))) {
      stop("Paper results must be unrenamed DStressR result tables.", call. = FALSE)
    }
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(schema_version = 1L, results = results, assays = assays,
               parameters = parameters, diagnostics = diagnostics,
               annotations = annotations, provenance = provenance), path)
  invisible(path)
}

read_paper_artifact <- function(path) {
  if (!file.exists(path)) {
    stop("Missing paper artifact: ", path, ". Run its analysis/diagnostic generator first.",
         call. = FALSE)
  }
  artifact <- read_analysis_rds(path)
  required <- c("results", "assays", "parameters", "diagnostics", "annotations", "provenance")
  if (!identical(artifact$schema_version, 1L) || !all(required %in% names(artifact))) {
    stop("Unsupported paper artifact: ", path, call. = FALSE)
  }
  artifact
}
