#' Simulate a chemical-genomics screen
#'
#' @param n_reporters,n_perturbations Dimensions excluding DMSO.
#' @param n_replicates Number of technical replicates.
#' @param sigma Observation noise standard deviation.
#' @param seed Optional random seed.
#' @return A data frame suitable for [prepare_assay()].
#' @export
simulate_screen <- function(n_reporters = 12, n_perturbations = 24, n_replicates = 2,
                            sigma = 0.15, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  reporters <- paste0("P", seq_len(n_reporters))
  perturbations <- c("DMSO", paste0("C", seq_len(n_perturbations)))
  grid <- expand.grid(
    reporter = reporters,
    perturbation = perturbations,
    replicate = paste0("r", seq_len(n_replicates)),
    stringsAsFactors = FALSE
  )
  grid$batch <- paste0("b", 1 + (as.integer(factor(grid$replicate)) %% 2))
  promoter_baseline <- stats::rnorm(n_reporters, 10, 0.7)
  names(promoter_baseline) <- reporters
  compound_global <- stats::rnorm(length(perturbations), 0, 0.25)
  names(compound_global) <- perturbations
  compound_global["DMSO"] <- 0
  specific <- matrix(0, nrow = n_reporters, ncol = length(perturbations),
                     dimnames = list(reporters, perturbations))
  specific[cbind(sample(reporters, 8, replace = TRUE),
                 sample(setdiff(perturbations, "DMSO"), 8, replace = TRUE))] <-
    stats::rnorm(8, 1.2, 0.25) * sample(c(-1, 1), 8, replace = TRUE)
  growth <- stats::rlnorm(nrow(grid), log(0.4), 0.2)
  eta <- promoter_baseline[grid$reporter] +
    log2(growth) +
    compound_global[grid$perturbation] +
    specific[cbind(grid$reporter, grid$perturbation)] +
    ifelse(grid$batch == "b2", 0.1, 0) +
    stats::rnorm(nrow(grid), 0, sigma)
  grid$od_16h.measured <- growth
  grid$LUX.AUC_16 <- 2^eta
  grid$truth_specific <- specific[cbind(grid$reporter, grid$perturbation)]
  grid$truth_global <- compound_global[grid$perturbation]
  grid
}
