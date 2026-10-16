# Gaussian null check with independently estimated reporter residual variances,
# shared reference means, and repeated use of each residual variance across contrasts.
source(file.path("analysis", "_helpers.R"))
out_dir <- analysis_output_dir("inference_audit_2026-10-16")
load_sources <- function(path) {
  e <- new.env(parent = .GlobalEnv)
  for (f in list.files(path, pattern = "[.]R$", full.names = TRUE)) sys.source(f, e)
  e
}
old <- load_sources(file.path(out_dir, "baseline_R"))
new <- load_sources(analysis_path("R"))
n_sim <- as.integer(Sys.getenv("DSTRESSR_NULL_SIMULATIONS", "1000"))
stopifnot(n_sim > 1L)
set.seed(20261016)
designs <- expand.grid(size = c("large_df", "small_df"),
                       noise = c("equal", "heterogeneous"), stringsAsFactors = FALSE)
rows <- list()
for (scenario in seq_len(nrow(designs))) {
  design <- designs[scenario, ]
  n <- 7L
  p <- if (design$size == "large_df") 94L else 6L
  n0 <- if (design$size == "large_df") 16L else 3L
  for (iteration in seq_len(n_sim)) {
    counts <- matrix(sample(if (design$size == "large_df") 6:8 else 2:3,
                             n * p, replace = TRUE), n, p)
    df <- n0 + rowSums(counts) - p - 1
    sigma2 <- if (design$noise == "equal") rep(1, n) else exp(rnorm(n, 0, 0.8))
    control_mean <- rnorm(n, sd = sqrt(sigma2 / n0))
    total <- matrix(rnorm(n * p), n, p) * sqrt(sigma2 / counts) - control_mean
    sample_sigma2 <- sigma2 * rchisq(n, df) / df
    v <- sample_sigma2 * (1 / counts + 1 / n0)
    delta <- total - rep(colMeans(total), each = n)
    centered_variance <- function(v) {
      ((n - 1) / n)^2 * v + (rep(colSums(v), each = n) - v) / n^2
    }
    raw_specific <- centered_variance(v)
    ws_df <- as.vector(vapply(seq_len(p), function(j) new$centered_contrast_df(v[, j], df), numeric(n)))
    focal_df <- rep(df, p)
    old_total <- old$eb_moderate_se(sqrt(as.vector(v)), focal_df)
    old_specific <- old$eb_moderate_se(sqrt(as.vector(centered_variance(matrix(old_total$se^2, n, p)))), focal_df)
    proposed <- new$eb_moderate_se(sqrt(as.vector(raw_specific)), ws_df)
    stats <- list(
      prechange = list(se = old_specific$se, df = old_specific$df),
      ordinary_focal = list(se = sqrt(as.vector(raw_specific)), df = focal_df),
      ordinary_WS = list(se = sqrt(as.vector(raw_specific)), df = ws_df),
      single_WS = list(se = proposed$se, df = proposed$df)
    )
    for (method in names(stats)) {
      s <- stats[[method]]
      pval <- 2 * pt(abs(as.vector(delta) / s$se), s$df, lower.tail = FALSE)
      by_reporter <- vapply(seq_len(n), function(a) any(p.adjust(pval[seq(a, n * p, n)], "BH") < 0.05), logical(1))
      rows[[length(rows) + 1L]] <- data.frame(
        size = design$size, noise = design$noise, iteration = iteration, method = method,
        raw_rejection_fraction = mean(pval < 0.05),
        reporter_family_rejection_fraction = mean(by_reporter),
        global_family_rejected = as.numeric(any(p.adjust(pval, "BH") < 0.05))
      )
    }
  }
  message("Completed null scenario ", scenario, "/", nrow(designs))
}
all <- do.call(rbind, rows)
groups <- split(all, interaction(all$size, all$noise, all$method, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(x) {
  out <- x[1, c("size", "noise", "method")]
  out$n_simulations <- nrow(x)
  for (field in c("raw_rejection_fraction", "reporter_family_rejection_fraction", "global_family_rejected")) {
    out[[field]] <- mean(x[[field]])
    out[[paste0(field, "_mc_se")]] <- sd(x[[field]]) / sqrt(nrow(x))
  }
  out
}))
utils::write.table(summary, file.path(out_dir, "null_calibration_summary.tsv"),
                    sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(all, file.path(out_dir, "null_calibration_repetitions.rds"))
print(summary, row.names = FALSE)
