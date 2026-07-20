analysis_project_root <- function() {
  env_root <- Sys.getenv("DSTRESSR_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }

  current <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(file.path(current, "R")) &&
        dir.exists(file.path(current, "analysis"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find the DStressR project root.", call. = FALSE)
    }
    current <- parent
  }
}

analysis_path <- function(...) {
  file.path(analysis_project_root(), ...)
}

analysis_output_dir <- function(...) {
  path <- analysis_path("analysis", "outputs", ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

analysis_data_root <- function() {
  env_root <- Sys.getenv("DSTRESSR_DATA_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }

  candidate <- analysis_path("workflow", "data")
  if (dir.exists(candidate)) {
    return(candidate)
  }

  local_candidates <- c(
    file.path(path.expand("~"), "Documents", "GitHub",
              "campylobacter_stressregnet", "workflow", "data"),
    file.path(dirname(analysis_project_root()),
              "campylobacter_stressregnet", "workflow", "data")
  )
  for (candidate in local_candidates) {
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    "External workflow data not found. Set DSTRESSR_DATA_ROOT to the ",
    "campylobacter_stressregnet/workflow/data directory.",
    call. = FALSE
  )
}

load_destress_package <- function() {
  project <- analysis_project_root()
  r_files <- list.files(file.path(project, "R"), pattern = "[.]R$", full.names = TRUE)
  for (path in r_files) {
    sys.source(path, envir = .GlobalEnv)
  }
  invisible(TRUE)
}

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

libmap_path <- function() {
  file.path(analysis_data_root(), "00-import", "Campylobacter", "LibMap.txt")
}

package_results_dir <- function(...) {
  analysis_path("analysis", "outputs", "package_results", ...)
}

comparison_results_dir <- function(...) {
  analysis_output_dir("comparisons", ...)
}

package_pair_result_path <- function(method) {
  package_results_dir(paste0(method, "_pair_results.tsv"))
}

read_package_pair_results <- function(method, path = package_pair_result_path(method)) {
  if (!file.exists(path)) {
    stop(
      "Missing DStressR package output for method `", method, "`: ", path,
      "\nExpected a TSV with columns: promoter, compound, effect, pvalue, ",
      "padj_global, padj_by_promoter.",
      "\nGenerate/export this table from the package first; analysis scripts ",
      "must not compute estimators or p-values.",
      call. = FALSE
    )
  }

  tab <- read_tsv_base(path)
  required <- c("promoter", "compound", "effect", "pvalue", "padj_global", "padj_by_promoter")
  missing <- setdiff(required, names(tab))
  if (length(missing)) {
    stop(
      "DStressR package output for method `", method, "` is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  tab <- tab[, required, drop = FALSE]
  tab$promoter <- as.character(tab$promoter)
  tab$compound <- as.character(tab$compound)
  tab$effect <- as.numeric(tab$effect)
  tab$pvalue <- as.numeric(tab$pvalue)
  tab$padj_global <- as.numeric(tab$padj_global)
  tab$padj_by_promoter <- as.numeric(tab$padj_by_promoter)
  tab$pair_id <- paste(tab$promoter, tab$compound, sep = "__")

  if (anyDuplicated(tab$pair_id)) {
    duplicated_pairs <- unique(tab$pair_id[duplicated(tab$pair_id)])
    stop(
      "DStressR package output for method `", method,
      "` is not one row per promoter-compound pair. First duplicated keys: ",
      paste(head(duplicated_pairs), collapse = ", "),
      call. = FALSE
    )
  }

  names(tab)[names(tab) %in% c("effect", "pvalue", "padj_global", "padj_by_promoter")] <-
    paste(method, c("effect", "pvalue", "padj_global", "padj_by_promoter"), sep = "_")
  tab
}

merge_package_pair_results <- function(methods) {
  if (length(methods) < 2) {
    stop("At least two methods are required for a comparison.", call. = FALSE)
  }

  tabs <- lapply(methods, read_package_pair_results)
  merged <- Reduce(
    function(x, y) merge(x, y, by = c("promoter", "compound", "pair_id"), all = FALSE, sort = FALSE),
    tabs
  )

  if (nrow(merged) == 0) {
    stop("No common promoter-compound pairs across package outputs.", call. = FALSE)
  }

  merged
}

safe_neglog10 <- function(x) {
  -log10(pmax(x, .Machine$double.xmin))
}

method_hit_column <- function(method) {
  paste0(method, "_hit")
}

comparison_adjustment <- function() {
  adjustment <- Sys.getenv("DSTRESSR_COMPARISON_ADJUSTMENT", unset = "global")
  adjustment <- gsub("-", "_", tolower(adjustment), fixed = TRUE)
  if (adjustment %in% c("promoter", "within_promoter")) {
    adjustment <- "by_promoter"
  }
  if (!adjustment %in% c("global", "by_promoter")) {
    stop(
      "DSTRESSR_COMPARISON_ADJUSTMENT must be `global` or `by_promoter`.",
      call. = FALSE
    )
  }
  adjustment
}

padj_column <- function(method, adjustment = comparison_adjustment()) {
  paste0(method, "_padj_", adjustment)
}

add_hit_columns <- function(tab, methods, fdr = 0.05, adjustment = comparison_adjustment()) {
  for (method in methods) {
    padj_col <- padj_column(method, adjustment)
    tab[[method_hit_column(method)]] <- is.finite(tab[[padj_col]]) & tab[[padj_col]] < fdr
  }
  tab
}

method_label <- function(method) {
  labels <- c(
    median_polish = "Median-polish max-p model",
    destress_standard = "DStressR ordinary model",
    destress_moderated = "DStressR default (moderated)",
    empty_vector = "Empty-vector control"
  )
  ifelse(method %in% names(labels), labels[[method]], method)
}
