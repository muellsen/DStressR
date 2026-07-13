make_key <- function(df, cols) {
  if (length(cols) == 0) {
    rep("all", nrow(df))
  } else {
    do.call(paste, c(df[, cols, drop = FALSE], sep = "\r"))
  }
}

permutation_null_effects <- function(null_values, center, n, B) {
  null_values <- null_values[is.finite(null_values)]
  if (length(null_values) < n || B <= 0) {
    return(NA_real_)
  }
  perm_means <- replicate(B, mean(sample(null_values, size = n, replace = FALSE)))
  perm_means - center
}

permutation_tail_pvalue <- function(effect, perm_effects, B, alternative) {
  if (!is.finite(effect) || all(!is.finite(perm_effects)) || B <= 0) {
    return(NA_real_)
  }
  if (alternative == "two.sided") {
    observed <- abs(effect)
    perm_stat <- abs(perm_effects)
    (1 + sum(perm_stat >= observed, na.rm = TRUE)) / (B + 1)
  } else if (alternative == "greater") {
    (1 + sum(perm_effects >= effect, na.rm = TRUE)) / (B + 1)
  } else {
    (1 + sum(perm_effects <= effect, na.rm = TRUE)) / (B + 1)
  }
}

#' Empirical p-values from replicate-averaged compound effects and DMSO nulls
#'
#' Computes empirical p-values by first averaging replicate-level values for
#' each promoter-compound-stratum combination, then comparing each averaged
#' non-control compound value to the corresponding distribution of averaged DMSO
#' control values from the same promoter and stratum, such as library plate.
#' Optionally, it also computes Monte Carlo permutation p-values by repeatedly
#' drawing replicate-sized DMSO sets from the same matched null stratum.
#'
#' @param table Replicate-level data frame.
#' @param value Numeric value column to test, for example a DStressR adjusted
#'   effect such as `destress_eb_effect_centered`.
#' @param promoter,compound Columns identifying promoters and compounds.
#' @param control Character vector of control compound IDs, usually DMSO wells.
#' @param replicate Optional replicate column. The function does not require
#'   replicate labels for averaging, but the argument documents the intended
#'   replicate-level input and is checked when supplied.
#' @param strata Optional columns defining matched null strata. Use
#'   `strata = "libplate"` to compare compounds only to DMSO wells from the
#'   same library plate.
#' @param min_replicates Minimum finite replicate values required for an
#'   averaged compound or DMSO control value.
#' @param min_null Minimum number of averaged DMSO controls required to compute
#'   an empirical p-value.
#' @param permutation If `TRUE`, compute an additional permutation p-value by
#'   sampling replicate-sized DMSO sets within each matched null stratum.
#' @param B Number of permutation draws when `permutation = TRUE`.
#' @param seed Optional random seed for reproducible permutation p-values.
#' @param alternative One of `"two.sided"`, `"greater"`, or `"less"`.
#' @param padj_method Multiple-testing correction method passed to
#'   [stats::p.adjust()], applied within promoter.
#' @return A data frame with one row per non-control promoter-compound-stratum
#'   average, including empirical p-values and promoter-wise adjusted p-values.
#' @export
empirical_replicate_pvalues <- function(table,
                                        value,
                                        promoter = "promoter",
                                        compound = "compound",
                                        control,
                                        replicate = NULL,
                                        strata = NULL,
                                        min_replicates = 2,
                                        min_null = 5,
                                        permutation = FALSE,
                                        B = 1000,
                                        seed = NULL,
                                        alternative = c("two.sided", "greater", "less"),
                                        padj_method = "BH") {
  stopifnot(is.data.frame(table))
  alternative <- match.arg(alternative)
  required <- c(value, promoter, compound, strata)
  if (!is.null(replicate)) {
    required <- c(required, replicate)
  }
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (missing(control)) {
    stop("Provide `control`, the control compound IDs used as empirical nulls.", call. = FALSE)
  }

  d <- table
  d$.empirical_value <- as.numeric(d[[value]])
  d <- d[is.finite(d$.empirical_value), , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No finite values available for empirical p-value computation.", call. = FALSE)
  }
  d[[promoter]] <- as.character(d[[promoter]])
  d[[compound]] <- as.character(d[[compound]])
  for (col in strata) {
    d[[col]] <- as.character(d[[col]])
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  average_cols <- unique(c(promoter, strata, compound))
  avg <- stats::aggregate(
    .empirical_value ~ .,
    d[, c(average_cols, ".empirical_value"), drop = FALSE],
    mean,
    na.rm = TRUE
  )
  names(avg)[names(avg) == ".empirical_value"] <- "mean_value"
  n_df <- stats::aggregate(
    .empirical_value ~ .,
    d[, c(average_cols, ".empirical_value"), drop = FALSE],
    function(x) sum(is.finite(x))
  )
  names(n_df)[names(n_df) == ".empirical_value"] <- "n_replicates"
  avg <- merge(avg, n_df, by = average_cols, all.x = TRUE, sort = FALSE)
  avg <- avg[avg$n_replicates >= min_replicates, , drop = FALSE]
  avg$is_control <- avg[[compound]] %in% control

  control_avg <- avg[avg$is_control, , drop = FALSE]
  test_avg <- avg[!avg$is_control, , drop = FALSE]
  if (nrow(control_avg) == 0) {
    stop("No averaged control rows available. Check `control` and `min_replicates`.", call. = FALSE)
  }
  if (nrow(test_avg) == 0) {
    stop("No averaged non-control rows available after filtering.", call. = FALSE)
  }

  null_cols <- unique(c(promoter, strata))
  d$.null_key <- make_key(d, null_cols)
  control_avg$.null_key <- make_key(control_avg, null_cols)
  test_avg$.null_key <- make_key(test_avg, null_cols)
  null_split <- split(control_avg$mean_value, control_avg$.null_key)
  raw_control_split <- split(
    d$.empirical_value[d[[compound]] %in% control],
    d$.null_key[d[[compound]] %in% control]
  )

  test_avg$null_n <- NA_integer_
  test_avg$null_center <- NA_real_
  test_avg$empirical_effect <- NA_real_
  test_avg$empirical_pvalue <- NA_real_
  test_avg$permutation_pvalue <- NA_real_
  test_avg$permutation_B <- if (isTRUE(permutation)) B else NA_integer_
  test_avg$permutation_min_pvalue <- if (isTRUE(permutation) && B > 0) 1 / (B + 1) else NA_real_
  permutation_cache <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(test_avg))) {
    null_values <- null_split[[test_avg$.null_key[i]]]
    null_values <- null_values[is.finite(null_values)]
    if (length(null_values) < min_null) {
      next
    }
    center <- stats::median(null_values, na.rm = TRUE)
    effect <- test_avg$mean_value[i] - center
    if (alternative == "two.sided") {
      score <- abs(effect)
      null_score <- abs(null_values - center)
      p <- (1 + sum(null_score >= score, na.rm = TRUE)) / (length(null_values) + 1)
    } else if (alternative == "greater") {
      p <- (1 + sum((null_values - center) >= effect, na.rm = TRUE)) / (length(null_values) + 1)
    } else {
      p <- (1 + sum((null_values - center) <= effect, na.rm = TRUE)) / (length(null_values) + 1)
    }
    test_avg$null_n[i] <- length(null_values)
    test_avg$null_center[i] <- center
    test_avg$empirical_effect[i] <- effect
    test_avg$empirical_pvalue[i] <- p
    if (isTRUE(permutation)) {
      cache_key <- paste(test_avg$.null_key[i], test_avg$n_replicates[i], sep = "\r")
      if (!exists(cache_key, envir = permutation_cache, inherits = FALSE)) {
        raw_null <- raw_control_split[[test_avg$.null_key[i]]]
        assign(
          cache_key,
          permutation_null_effects(
            null_values = raw_null,
            center = center,
            n = test_avg$n_replicates[i],
            B = B
          ),
          envir = permutation_cache
        )
      }
      perm_effects <- get(cache_key, envir = permutation_cache, inherits = FALSE)
      test_avg$permutation_pvalue[i] <- permutation_tail_pvalue(
        effect = effect,
        perm_effects = perm_effects,
        B = B,
        alternative = alternative
      )
    }
  }

  test_avg$empirical_padj_by_promoter <- NA_real_
  test_avg$permutation_padj_by_promoter <- NA_real_
  for (idx in split(seq_len(nrow(test_avg)), test_avg[[promoter]])) {
    ok <- idx[is.finite(test_avg$empirical_pvalue[idx])]
    test_avg$empirical_padj_by_promoter[ok] <- stats::p.adjust(
      test_avg$empirical_pvalue[ok],
      method = padj_method
    )
    ok_perm <- idx[is.finite(test_avg$permutation_pvalue[idx])]
    test_avg$permutation_padj_by_promoter[ok_perm] <- stats::p.adjust(
      test_avg$permutation_pvalue[ok_perm],
      method = padj_method
    )
  }
  test_avg$.null_key <- NULL
  test_avg$is_control <- NULL
  test_avg
}
