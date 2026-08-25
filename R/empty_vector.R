#' Reproduce Empty Vector Control normalization
#'
#' This function implements the Salmonella StressRegNet workflow in which
#' reporter activity is normalized against an Empty Vector Control (EVC) reporter
#' measured for the same perturbation/library well. It starts from a long expression
#' table, subtracts a perturbation-specific EVC average from each
#' reporter-replicate value, estimates reporter-replicate DMSO null
#' distributions, and applies the original conservative replicate aggregation.
#'
#' @param data Long expression table with one row per
#'   reporter-perturbation-replicate observation.
#' @param reporter,perturbation,replicate Column names identifying reporter,
#'   perturbation/library well, and replicate.
#' @param response Column containing the expression value to normalize. For the
#'   Salmonella workflow this is `log2.lux.normed.centered`.
#' @param empty_vector_reporter Reporter/control strain used as the Empty Vector
#'   reference. The original Salmonella workflow uses `PEVC3`.
#' @param control Character vector of perturbation/library-well IDs used as DMSO
#'   controls for the null distribution.
#' @param exclude Character vector of perturbation/library-well IDs removed before
#'   normalization and hit calling, for example noisy DMSO wells.
#' @param remove_reporters Reporters removed before normalization, for example
#'   failed reporter strains.
#' @param fdr FDR threshold used to assign the `hit` class in the pair-level
#'   table.
#' @param require_complete_empty_vector If `TRUE`, require all EVC replicate
#'   values for a perturbation to be finite before computing the EVC average. This
#'   matches the original workflow's effective behavior with two PEVC3
#'   replicates.
#' @return A list of class `destress_empty_vector` with `replicate_results`,
#'   `pair_results`, `empty_vector_reference`, `control`, and `exclude`
#'   components.
#' @export
fit_empty_vector_control <- function(data,
                                     reporter = "reporter",
                                     perturbation = "srn_code",
                                     replicate = "replicate",
                                     response = "log2.lux.normed.centered",
                                     empty_vector_reporter = "PEVC3",
                                     control,
                                     exclude = character(),
                                     remove_reporters = character(),
                                     fdr = 0.05,
                                     require_complete_empty_vector = TRUE) {
  stopifnot(is.data.frame(data))
  if (missing(control) || length(control) == 0) {
    stop("`control` must contain one or more DMSO/control perturbation IDs.", call. = FALSE)
  }
  required <- c(reporter, perturbation, replicate, response)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- data
  d$.reporter <- as.character(d[[reporter]])
  d$.perturbation <- as.character(d[[perturbation]])
  d$.replicate <- as.character(d[[replicate]])
  d$.response <- as.numeric(d[[response]])
  d <- d[is.finite(d$.response), , drop = FALSE]
  d <- d[!(d$.perturbation %in% exclude), , drop = FALSE]
  d <- d[!(d$.reporter %in% remove_reporters), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite rows remain after applying exclusions.", call. = FALSE)
  }

  ev <- d[d$.reporter == empty_vector_reporter, , drop = FALSE]
  if (nrow(ev) == 0) {
    stop("No rows found for `empty_vector_reporter = \"", empty_vector_reporter, "\"`.", call. = FALSE)
  }
  ev_counts <- stats::aggregate(
    .response ~ .perturbation,
    ev,
    function(x) c(
      empty_vector_mean = mean(x, na.rm = TRUE),
      empty_vector_n = sum(is.finite(x)),
      empty_vector_missing = sum(!is.finite(x))
    )
  )
  ev_counts <- do.call(data.frame, ev_counts)
  names(ev_counts) <- c(
    "srn_code",
    "empty_vector_mean",
    "empty_vector_n",
    "empty_vector_missing"
  )
  if (isTRUE(require_complete_empty_vector)) {
    expected_n <- max(ev_counts$empty_vector_n, na.rm = TRUE)
    ev_counts$empty_vector_mean[ev_counts$empty_vector_n < expected_n] <- NA_real_
  }

  d <- merge(
    d,
    ev_counts,
    by.x = ".perturbation",
    by.y = "srn_code",
    all.x = TRUE,
    sort = FALSE
  )
  d$log.evcfc <- d$.response - d$empty_vector_mean
  d <- d[is.finite(d$log.evcfc) & d$.reporter != empty_vector_reporter, , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite Empty Vector-normalized rows remain.", call. = FALSE)
  }
  if (!any(d$.perturbation %in% control)) {
    stop("No control rows found after Empty Vector normalization.", call. = FALSE)
  }

  d$promoter_replicate <- paste(d$.reporter, d$.replicate, sep = "_")
  dmso <- d[d$.perturbation %in% control, , drop = FALSE]
  dmso_params <- stats::aggregate(
    log.evcfc ~ promoter_replicate,
    dmso,
    function(x) c(
      dmso.mean = mean(x, na.rm = TRUE),
      dmso.stdv = stats::sd(x, na.rm = TRUE),
      n = sum(is.finite(x))
    )
  )
  dmso_params <- do.call(data.frame, dmso_params)
  names(dmso_params) <- c("promoter_replicate", "dmso.mean", "dmso.stdv", "n")

  replicate_results <- merge(d, dmso_params, by = "promoter_replicate", all.x = TRUE, sort = FALSE)
  replicate_results$zscore <- (
    replicate_results$log.evcfc - replicate_results$dmso.mean
  ) / replicate_results$dmso.stdv
  replicate_results$pvalue <- 2 * stats::pnorm(
    q = abs(replicate_results$zscore),
    mean = 0,
    sd = 1,
    lower.tail = FALSE
  )
  replicate_results <- replicate_results[, c(
    "promoter_replicate",
    ".reporter",
    ".replicate",
    ".perturbation",
    "log.evcfc",
    "empty_vector_mean",
    "empty_vector_n",
    "dmso.mean",
    "dmso.stdv",
    "n",
    "zscore",
    "pvalue"
  )]
  names(replicate_results)[names(replicate_results) == ".reporter"] <- "reporter"
  names(replicate_results)[names(replicate_results) == ".replicate"] <- "replicate"
  names(replicate_results)[names(replicate_results) == ".perturbation"] <- "srn_code"

  pair_results <- replicate_results[is.finite(replicate_results$pvalue), , drop = FALSE]
  pair_results <- pair_results[order(
    pair_results$reporter,
    pair_results$srn_code,
    -pair_results$pvalue
  ), , drop = FALSE]
  pair_results <- pair_results[!duplicated(paste(pair_results$reporter, pair_results$srn_code)), , drop = FALSE]
  pair_results$pvalue.adj <- NA_real_
  split_idx <- split(seq_len(nrow(pair_results)), pair_results$reporter)
  for (idx in split_idx) {
    pair_results$pvalue.adj[idx] <- stats::p.adjust(pair_results$pvalue[idx], method = "BH")
  }
  pair_results$hit <- "Not DE"
  up <- is.finite(pair_results$pvalue.adj) &
    pair_results$pvalue.adj < fdr &
    pair_results$log.evcfc > 0
  down <- is.finite(pair_results$pvalue.adj) &
    pair_results$pvalue.adj < fdr &
    pair_results$log.evcfc < 0
  pair_results$hit[up] <- "Upregulated"
  pair_results$hit[down] <- "Downregulated"

  structure(
    list(
      replicate_results = replicate_results,
      pair_results = pair_results,
      empty_vector_reference = ev_counts,
      control = control,
      exclude = exclude,
      empty_vector_reporter = empty_vector_reporter
    ),
    class = "destress_empty_vector"
  )
}
