source(file.path("analysis", "_helpers.R"))
source(file.path("analysis", "_paper_artifacts.R"))
load_destress_package()

load_binsfeld_paper_data()
out_dir <- analysis_output_dir("binsfeld_modeling_steps")
evc_method_dir <- analysis_output_dir("binsfeld_evc_calibrated")

wt_auc <- binsfeld_reporter_auc[
  binsfeld_reporter_auc$strain == "WT" &
    binsfeld_reporter_auc$removed == "No",
]

prepare_binsfeld <- function(growth_exponent) {
  prepare_assay(
    wt_auc,
    reporter = "promoter",
    perturbation = "compound",
    control = "Water",
    lux = "lux_auc",
    growth = "od_auc",
    growth_exponent = growth_exponent,
    batch = "dose_level",
    replicate = "replicate",
    growth_covariates = "replicate",
    numeric_covariates = "dose_level"
  )
}

modeled <- prepare_binsfeld("estimate")
alpha1 <- prepare_binsfeld(1)
evc_huber <- prepare_assay(
  wt_auc,
  reporter = "promoter",
  perturbation = "compound",
  control = "Water",
  lux = "lux_auc",
  growth = "od_auc",
  growth_exponent = "estimate",
  batch = "dose_level",
  replicate = "replicate",
  growth_covariates = "replicate",
  numeric_covariates = "dose_level",
  background_reporter = "EVC",
  background_method = "huber",
  background_by = c("compound", "dose_level", "replicate")
)

growth_parameters <- attr(modeled, "destress")$growth_exponent_fit
utils::write.table(
  growth_parameters,
  file.path(out_dir, "binsfeld_growth_parameter_estimates.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

background_calibration <- attr(evc_huber, "destress")$background_fit
if (!is.null(background_calibration)) {
  utils::write.table(
    background_calibration,
    file.path(out_dir, "binsfeld_evc_background_calibration.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

response_matrix <- function(assay) {
  tab <- assay[, c("promoter", "compound", ".response")]
  means <- stats::aggregate(.response ~ promoter + compound, tab, mean, na.rm = TRUE)
  control <- means[means$compound == "Water", c("promoter", ".response")]
  names(control)[2] <- "control_response"
  means <- merge(means, control, by = "promoter", all.x = TRUE, sort = FALSE)
  means$delta_response <- means$.response - means$control_response
  means <- means[means$compound != "Water" & means$promoter != "EVC", ]
  means
}

estimated_exponent_matrix <- response_matrix(modeled)
alpha1_matrix <- response_matrix(alpha1)
evc_huber_matrix <- response_matrix(evc_huber)
names(estimated_exponent_matrix)[names(estimated_exponent_matrix) == "delta_response"] <- "estimated_exponent_response"
names(alpha1_matrix)[names(alpha1_matrix) == "delta_response"] <- "alpha1_response"
names(evc_huber_matrix)[names(evc_huber_matrix) == "delta_response"] <- "evc_huber_response"

score_keys <- unique(wt_auc[, c("well", "drug", "promoter", "replicate", "compound")])
binsfeld_scores <- binsfeld_reporter_scores[
  binsfeld_reporter_scores$strain == "WT" &
    binsfeld_reporter_scores$statistic == "Scores",
  ,
  drop = FALSE
]
binsfeld_scores <- merge(
  binsfeld_scores[, c("well", "drug", "promoter", "replicate", "value")],
  score_keys,
  by = c("well", "drug", "promoter", "replicate"),
  all = FALSE,
  sort = FALSE
)
binsfeld_adjusted_response_matrix <- stats::aggregate(
  value ~ promoter + compound,
  binsfeld_scores[
    binsfeld_scores$compound != "Water" &
      binsfeld_scores$promoter != "EVC" &
      is.finite(binsfeld_scores$value),
    ,
    drop = FALSE
  ],
  mean
)
names(binsfeld_adjusted_response_matrix)[3] <- "binsfeld_adjusted_response"

matched <- merge(
  estimated_exponent_matrix[, c("promoter", "compound", "estimated_exponent_response")],
  alpha1_matrix[, c("promoter", "compound", "alpha1_response")],
  by = c("promoter", "compound"),
  all = FALSE,
  sort = FALSE
)
matched$difference <- matched$estimated_exponent_response - matched$alpha1_response

response_construction <- Reduce(
  function(x, y) merge(x, y, by = c("promoter", "compound"), all = FALSE, sort = FALSE),
  list(
    alpha1_matrix[, c("promoter", "compound", "alpha1_response")],
    estimated_exponent_matrix[, c("promoter", "compound", "estimated_exponent_response")],
    evc_huber_matrix[, c("promoter", "compound", "evc_huber_response")],
    binsfeld_adjusted_response_matrix[, c("promoter", "compound", "binsfeld_adjusted_response")]
  )
)
response_construction$estimated_exponent_minus_alpha1 <- response_construction$estimated_exponent_response -
  response_construction$alpha1_response
response_construction$evc_huber_minus_alpha1 <- response_construction$evc_huber_response -
  response_construction$alpha1_response
response_scatter_data <- response_construction

union_file <- file.path(evc_method_dir, "evc_huber_comparison_to_binsfeld_and_default.tsv")
if (file.exists(union_file)) {
  union_pairs <- read.delim(union_file, check.names = FALSE, stringsAsFactors = FALSE)
  union_pairs <- union_pairs[
    union_pairs$binsfeld_hit | union_pairs$estimated_exponent_hit | union_pairs$evc_huber_hit,
    ,
    drop = FALSE
  ]
  keep_compounds <- sort(unique(union_pairs$compound))
  matched <- matched[matched$compound %in% keep_compounds, ]
  response_construction <- response_construction[response_construction$compound %in% keep_compounds, ]
}

matrix_wide <- stats::reshape(
  matched[, c("promoter", "compound", "estimated_exponent_response")],
  idvar = "promoter",
  timevar = "compound",
  direction = "wide"
)
compound_cols <- setdiff(names(matrix_wide), "promoter")
compound_order <- sub("^estimated_exponent_response[.]", "", compound_cols)
if (length(compound_cols) > 2) {
  mat <- as.matrix(matrix_wide[, compound_cols, drop = FALSE])
  rownames(mat) <- matrix_wide$promoter
  mat[is.na(mat)] <- 0
  compound_order <- sub("^estimated_exponent_response[.]", "", colnames(mat)[stats::hclust(stats::dist(t(mat)))$order])
}
write_paper_artifact(file.path(out_dir, "binsfeld_response_paper.rds"),
  assays = list(modeled = modeled, alpha1 = alpha1, with_EVC = evc_huber),
  parameters = list(growth_parameters = growth_parameters, background_calibration = background_calibration),
  diagnostics = list(matched = matched, response_construction = response_construction,
                     response_scatter_data = response_scatter_data, compound_order = compound_order))
