#' Reproduce the legacy median-polish workflow
#'
#' This function implements the original median-polish hit-determination
#' workflow used for the Campylobacter promoter-library screen. It starts from a
#' long expression table, centers each promoter-library-plate-replicate group by
#' its DMSO wells, applies [stats::medpolish()] to the resulting
#' promoter-libplate-replicate by compound matrix, and computes z-test p-values
#' from the polished DMSO residual distribution.
#'
#' The promoter-compound hit table follows the original conservative replicate
#' aggregation: DMSO and excluded control wells are removed, the largest
#' replicate-level p-value is retained for each promoter-compound pair, and
#' p-values are BH-adjusted within promoter.
#'
#' @param data Long expression table with one row per
#'   promoter-compound-replicate observation.
#' @param promoter,compound,libplate,replicate Column names identifying the
#'   promoter, compound/library well, library plate, and replicate.
#' @param response Column containing the already growth-normalized log2
#'   response, for example `log2.auc.16hmeasured.normed`.
#' @param control Character vector of compound/library-well IDs used as DMSO
#'   controls.
#' @param exclude Character vector of compound/library-well IDs to remove before
#'   median polishing and hit calling, for example noisy DMSO wells.
#' @param fdr FDR threshold used to assign the `hit` class in the pair-level
#'   table.
#' @param normality If `TRUE`, test pre-polish DMSO-centered fold changes
#'   within each promoter-library-plate-replicate group.
#' @param normality_methods Character vector containing `"shapiro"` and/or
#'   `"lilliefors"`. The Lilliefors test requires the suggested `nortest`
#'   package.
#' @param maxiter,eps Passed to [stats::medpolish()].
#' @return A list of class `destress_median_polish` with `replicate_results`,
#'   `pair_results`, `polished_matrix`, `medpolish`, and optional
#'   `normality_results` components.
#' @export
fit_median_polish <- function(data,
                              promoter = "promoter",
                              compound = "srn_code",
                              libplate = "libplate",
                              replicate = "replicate",
                              response = "log2.auc.16hmeasured.normed",
                              control,
                              exclude = character(),
                              fdr = 0.05,
                              normality = FALSE,
                              normality_methods = c("shapiro", "lilliefors"),
                              maxiter = 1000,
                              eps = 1e-8) {
  stopifnot(is.data.frame(data))
  if (missing(control) || length(control) == 0) {
    stop("`control` must contain one or more DMSO/control compound IDs.", call. = FALSE)
  }
  required <- c(promoter, compound, libplate, replicate, response)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- data
  d$.promoter <- as.character(d[[promoter]])
  d$.compound <- as.character(d[[compound]])
  d$.libplate <- as.character(d[[libplate]])
  d$.replicate <- as.character(d[[replicate]])
  d$.response <- as.numeric(d[[response]])
  d$.group <- paste(d$.promoter, d$.libplate, d$.replicate, sep = "_")
  d <- d[is.finite(d$.response), , drop = FALSE]
  d <- d[!(d$.compound %in% exclude), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite rows remain after applying exclusions.", call. = FALSE)
  }
  if (!any(d$.compound %in% control)) {
    stop("No control rows found after applying exclusions.", call. = FALSE)
  }

  dmso <- d[d$.compound %in% control, , drop = FALSE]
  dmso_mean <- stats::aggregate(
    .response ~ .group,
    dmso,
    mean,
    na.rm = TRUE
  )
  names(dmso_mean)[2] <- ".dmso_mean"
  d <- merge(d, dmso_mean, by = ".group", all.x = TRUE, sort = FALSE)
  d$.log2FC <- d$.response - d$.dmso_mean
  d <- d[is.finite(d$.log2FC), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No rows have a finite DMSO-centered log2FC.", call. = FALSE)
  }
  normality_results <- NULL
  if (isTRUE(normality)) {
    normality_results <- dmso_normality_tests(
      d,
      control = control,
      methods = normality_methods
    )
  }

  key <- paste(d$.group, d$.compound, sep = "\r")
  if (anyDuplicated(key)) {
    stop(
      "Median-polish input must have at most one row per ",
      "promoter/libplate/replicate/compound combination.",
      call. = FALSE
    )
  }

  groups <- sort(unique(d$.group))
  compounds <- sort(unique(d$.compound))
  mat <- matrix(
    NA_real_,
    nrow = length(groups),
    ncol = length(compounds),
    dimnames = list(groups, compounds)
  )
  mat[cbind(match(d$.group, groups), match(d$.compound, compounds))] <- d$.log2FC

  polish <- stats::medpolish(mat, na.rm = TRUE, maxiter = maxiter, eps = eps, trace.iter = FALSE)
  polished <- polish$residuals
  long <- as.data.frame(as.table(polished), stringsAsFactors = FALSE)
  names(long) <- c("promoter_libplate_replicate", "srn_code", "log2FC.polished")
  long <- long[is.finite(long$log2FC.polished), , drop = FALSE]

  group_map <- unique(d[, c(".group", ".promoter", ".libplate", ".replicate"), drop = FALSE])
  names(group_map) <- c("promoter_libplate_replicate", "promoter", "libplate", "replicate")
  long <- merge(long, group_map, by = "promoter_libplate_replicate", all.x = TRUE, sort = FALSE)

  dmso_polished <- long[long$srn_code %in% control, , drop = FALSE]
  dmso_params <- stats::aggregate(
    log2FC.polished ~ promoter_libplate_replicate,
    dmso_polished,
    function(x) c(
      dmso.avg_dmsoFC = mean(x, na.rm = TRUE),
      dmso.stdv_dmsoFC = stats::sd(x, na.rm = TRUE),
      n = sum(is.finite(x))
    )
  )
  dmso_params <- do.call(data.frame, dmso_params)
  names(dmso_params) <- c(
    "promoter_libplate_replicate",
    "dmso.avg_dmsoFC",
    "dmso.stdv_dmsoFC",
    "n"
  )

  replicate_results <- merge(
    long,
    dmso_params,
    by = "promoter_libplate_replicate",
    all.x = TRUE,
    sort = FALSE
  )
  replicate_results$zscore <- (
    replicate_results$log2FC.polished - replicate_results$dmso.avg_dmsoFC
  ) / replicate_results$dmso.stdv_dmsoFC
  replicate_results$pvalue <- 2 * stats::pnorm(
    q = abs(replicate_results$zscore),
    mean = 0,
    sd = 1,
    lower.tail = FALSE
  )
  replicate_results <- replicate_results[, c(
    "promoter_libplate_replicate",
    "promoter",
    "libplate",
    "replicate",
    "srn_code",
    "log2FC.polished",
    "dmso.avg_dmsoFC",
    "dmso.stdv_dmsoFC",
    "n",
    "zscore",
    "pvalue"
  )]

  pair_results <- replicate_results[
    !(replicate_results$srn_code %in% c(control, exclude)) &
      is.finite(replicate_results$pvalue),
    ,
    drop = FALSE
  ]
  pair_results <- pair_results[order(
    pair_results$promoter,
    pair_results$srn_code,
    -pair_results$pvalue
  ), , drop = FALSE]
  pair_results <- pair_results[!duplicated(paste(pair_results$promoter, pair_results$srn_code)), , drop = FALSE]
  pair_results$pvalue.adj <- rep(NA_real_, nrow(pair_results))
  if (nrow(pair_results) > 0) {
    split_idx <- split(seq_len(nrow(pair_results)), pair_results$promoter)
    for (idx in split_idx) {
      pair_results$pvalue.adj[idx] <- stats::p.adjust(pair_results$pvalue[idx], method = "BH")
    }
  }
  pair_results$hit <- rep("Not DE", nrow(pair_results))
  up <- is.finite(pair_results$pvalue.adj) &
    pair_results$pvalue.adj < fdr &
    pair_results$log2FC.polished > 0
  down <- is.finite(pair_results$pvalue.adj) &
    pair_results$pvalue.adj < fdr &
    pair_results$log2FC.polished < 0
  pair_results$hit[up] <- "Upregulated"
  pair_results$hit[down] <- "Downregulated"

  structure(
    list(
      replicate_results = replicate_results,
      pair_results = pair_results,
      polished_matrix = polished,
      medpolish = polish,
      control = control,
      exclude = exclude,
      normality_results = normality_results
    ),
    class = "destress_median_polish"
  )
}

dmso_normality_tests <- function(data, control, methods = c("shapiro", "lilliefors")) {
  methods <- gsub("-", "_", tolower(methods), fixed = TRUE)
  methods[methods == "lillie"] <- "lilliefors"
  unknown <- setdiff(methods, c("shapiro", "lilliefors"))
  if (length(unknown) > 0) {
    stop("Unknown normality methods: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if ("lilliefors" %in% methods && !requireNamespace("nortest", quietly = TRUE)) {
    stop(
      "The Lilliefors normality test requires the `nortest` package. ",
      "Install it or use `normality_methods = \"shapiro\"`.",
      call. = FALSE
    )
  }

  dmso <- data[data$.compound %in% control, , drop = FALSE]
  groups <- split(dmso, dmso$.group)
  out <- do.call(rbind, lapply(groups, function(group_data) {
    x <- group_data$.log2FC[is.finite(group_data$.log2FC)]
    shapiro_p <- NA_real_
    lillie_p <- NA_real_
    if (length(x) >= 3 && length(x) <= 5000 && stats::sd(x) > 0) {
      if ("shapiro" %in% methods) {
        shapiro_p <- stats::shapiro.test(x)$p.value
      }
      if ("lilliefors" %in% methods) {
        lillie_p <- nortest::lillie.test(x)$p.value
      }
    }
    data.frame(
      promoter_libplate_replicate = group_data$.group[1],
      promoter = group_data$.promoter[1],
      libplate = group_data$.libplate[1],
      replicate = group_data$.replicate[1],
      n = length(x),
      shapiro.pval = shapiro_p,
      lillie.pval = lillie_p,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  if ("shapiro" %in% methods) {
    out$shapiro.pval.adj <- stats::p.adjust(out$shapiro.pval, method = "BH")
  }
  if ("lilliefors" %in% methods) {
    out$lillie.pval.adj <- stats::p.adjust(out$lillie.pval, method = "BH")
  }
  out
}
