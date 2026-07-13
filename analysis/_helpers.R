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
