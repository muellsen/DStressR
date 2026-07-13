#' Simulate a chemical-genomics screen
#'
#' @param n_promoters,n_compounds Dimensions excluding DMSO.
#' @param n_replicates Number of technical replicates.
#' @param sigma Observation noise standard deviation.
#' @param seed Optional random seed.
#' @return A data frame suitable for [prepare_assay()].
#' @export
simulate_screen <- function(n_promoters = 12, n_compounds = 24, n_replicates = 2,
                            sigma = 0.15, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  promoters <- paste0("P", seq_len(n_promoters))
  compounds <- c("DMSO", paste0("C", seq_len(n_compounds)))
  grid <- expand.grid(
    promoter = promoters,
    compound = compounds,
    replicate = paste0("r", seq_len(n_replicates)),
    stringsAsFactors = FALSE
  )
  grid$batch <- paste0("b", 1 + (as.integer(factor(grid$replicate)) %% 2))
  promoter_baseline <- stats::rnorm(n_promoters, 10, 0.7)
  names(promoter_baseline) <- promoters
  compound_global <- stats::rnorm(length(compounds), 0, 0.25)
  names(compound_global) <- compounds
  compound_global["DMSO"] <- 0
  specific <- matrix(0, nrow = n_promoters, ncol = length(compounds),
                     dimnames = list(promoters, compounds))
  specific[cbind(sample(promoters, 8, replace = TRUE),
                 sample(setdiff(compounds, "DMSO"), 8, replace = TRUE))] <-
    stats::rnorm(8, 1.2, 0.25) * sample(c(-1, 1), 8, replace = TRUE)
  growth <- stats::rlnorm(nrow(grid), log(0.4), 0.2)
  eta <- promoter_baseline[grid$promoter] +
    log2(growth) +
    compound_global[grid$compound] +
    specific[cbind(grid$promoter, grid$compound)] +
    ifelse(grid$batch == "b2", 0.1, 0) +
    stats::rnorm(nrow(grid), 0, sigma)
  grid$od_16h.measured <- growth
  grid$LUX.AUC_16 <- 2^eta
  grid$truth_specific <- specific[cbind(grid$promoter, grid$compound)]
  grid$truth_global <- compound_global[grid$compound]
  grid
}
