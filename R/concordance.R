# ===================================================================== #
#  Licensed as GPL-v3.0.                                                #
#                                                                       #
#  Developed as part of the AMRverse (https://github.com/AMRverse):     #
#  https://github.com/AMRverse/AMRgen                                   #
#                                                                       #
#  We created this package for both routine data analysis and academic  #
#  research and it was publicly released in the hope that it will be    #
#  useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                       #
#  This R package is free software; you can freely use and distribute   #
#  it for both personal and commercial purposes under the terms of the  #
#  GNU General Public License version 3.0 (GNU GPL-3), as published by  #
#  the Free Software Foundation.                                        #
# ===================================================================== #

#' Calculate Genotype-Phenotype Concordance
#'
#' Compares genotypic predictions (presence of resistance markers) to phenotypic
#' truth (resistant vs susceptible) using a binary matrix from [get_binary_matrix()].
#' Calculates standard classification metrics via yardstick and AMR-specific error
#' rates (VME and ME) per ISO 20776-2.
#'
#' @param binary_matrix A data frame output by [get_binary_matrix()], containing
#'   sample-level binary marker presence/absence alongside phenotype columns
#'   (`R`, `NWT`).
#' @param markers A character vector of marker column names to include in the
#'   genotypic prediction. Default `NULL` includes all marker columns.
#' @param exclude_markers A character vector of marker column names to exclude
#'   from the genotypic prediction. Applied after `markers` filtering.
#' @param ppv_threshold A numeric PPV threshold (0-1). Markers with solo PPV
#'   below this value are excluded. Requires `solo_ppv_results`.
#' @param solo_ppv_results Output of [solo_ppv_analysis()], used for PPV-based
#'   marker filtering when `ppv_threshold` is set.
#' @param truth A character string specifying the phenotypic truth column to use:
#'   `"R"` (resistant vs susceptible/intermediate, default) or `"NWT"`
#'   (non-wildtype vs wildtype).
#' @param prediction_rule A character string specifying the rule for generating
#'   genotypic predictions. Currently only `"any"` is supported: a sample is
#'   predicted positive if any included marker is present (value of 1).
#'
#' @details
#' The function identifies marker columns as all columns not in the reserved set
#' (`id`, `pheno`, `ecoff`, `R`, `I`, `NWT`, `mic`, `disk`). It then applies
#' filtering in order: custom `markers` list, `exclude_markers`, then
#' `ppv_threshold` filtering.
#'
#' Genotypic prediction is generated per sample: if `prediction_rule = "any"`,
#' then a sample is predicted positive (1) if any selected marker equals 1.
#'
#' Standard metrics (sensitivity, specificity, PPV, NPV, accuracy, kappa,
#' F-measure) are calculated using yardstick. AMR-specific error rates are
#' computed internally:
#' - **VME** (Very Major Error): FN / (TP + FN) = 1 - sensitivity. Proportion of
#'   truly resistant isolates missed by genotype.
#' - **ME** (Major Error): FP / (TN + FP) = 1 - specificity. Proportion of
#'   truly susceptible isolates incorrectly called resistant by genotype.
#'
#' @return An S3 object of class `"amr_concordance"`, a list containing:
#' - `conf_mat`: A yardstick confusion matrix object.
#' - `metrics`: A tibble with columns `.metric`, `.estimator`, `.estimate`
#'   containing sensitivity, specificity, ppv, npv, accuracy, kap, f_meas, VME,
#'   and ME.
#' - `data`: The input binary matrix with an added `geno_prediction` column.
#' - `markers_used`: Character vector of markers included in the prediction.
#' - `truth_col`: Which truth column was used (`"R"` or `"NWT"`).
#' - `n`: Total number of samples with non-missing truth values.
#'
#' @importFrom dplyr across any_of filter mutate select
#' @importFrom tibble tibble
#' @importFrom yardstick conf_mat sens spec ppv npv accuracy kap f_meas
#' @seealso [get_binary_matrix()], [solo_ppv_analysis()]
#' @export
#' @examples
#' \dontrun{
#' geno_table <- import_amrfp(ecoli_geno_raw, "Name")
#'
#' binary_matrix <- get_binary_matrix(
#'   geno_table = geno_table,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno_clsi"
#' )
#'
#' # Basic concordance using all markers
#' result <- concordance(binary_matrix)
#' result
#'
#' # Exclude specific markers
#' result <- concordance(binary_matrix, exclude_markers = c("qnrS1"))
#'
#' # Filter markers by solo PPV threshold
#' solo_ppv <- solo_ppv_analysis(binary_matrix = binary_matrix)
#' result <- concordance(
#'   binary_matrix,
#'   ppv_threshold = 0.5,
#'   solo_ppv_results = solo_ppv
#' )
#'
#' # Access components
#' result$conf_mat
#' result$metrics
#' result$markers_used
#' }
concordance <- function(binary_matrix,
                        markers = NULL,
                        exclude_markers = NULL,
                        ppv_threshold = NULL,
                        solo_ppv_results = NULL,
                        truth = "R",
                        prediction_rule = "any") {
  # --- input validation ---
  if (!is.data.frame(binary_matrix)) {
    stop("`binary_matrix` must be a data frame (output of get_binary_matrix()).")
  }

  truth <- match.arg(truth, choices = c("R", "NWT"))
  prediction_rule <- match.arg(prediction_rule, choices = c("any"))

  if (!(truth %in% colnames(binary_matrix))) {
    stop(paste0("Truth column '", truth, "' not found in binary_matrix."))
  }

  if (!is.null(ppv_threshold) && is.null(solo_ppv_results)) {
    stop("`solo_ppv_results` must be provided when `ppv_threshold` is set.")
  }

  # --- identify marker columns ---
  reserved_cols <- c("id", "pheno", "ecoff", "R", "I", "NWT", "mic", "disk")
  all_markers <- setdiff(colnames(binary_matrix), reserved_cols)

  if (length(all_markers) == 0) {
    stop("No marker columns found in binary_matrix.")
  }

  # --- apply marker filtering ---
  selected_markers <- all_markers

  # custom markers list
  if (!is.null(markers)) {
    unknown <- setdiff(markers, all_markers)
    if (length(unknown) > 0) {
      warning(paste("Markers not found in binary_matrix (ignored):",
                    paste(unknown, collapse = ", ")))
    }
    selected_markers <- intersect(markers, all_markers)
  }

  # exclude markers
  if (!is.null(exclude_markers)) {
    selected_markers <- setdiff(selected_markers, exclude_markers)
  }

  # PPV threshold filter
  if (!is.null(ppv_threshold)) {
    ppv_category <- if (truth == "R") "R" else "NWT"
    ppv_data <- solo_ppv_results$solo_stats
    passing_markers <- ppv_data %>%
      filter(category == ppv_category, ppv >= ppv_threshold) %>%
      pull(marker)
    selected_markers <- intersect(selected_markers, passing_markers)
  }

  if (length(selected_markers) == 0) {
    stop("No markers remaining after filtering. Adjust `markers`, `exclude_markers`, or `ppv_threshold`.")
  }

  # --- generate genotypic prediction ---
  df <- binary_matrix %>%
    filter(!is.na(get(truth)))

  if (nrow(df) == 0) {
    stop(paste0("No samples with non-NA values in truth column '", truth, "'."))
  }

  if (prediction_rule == "any") {
    df <- df %>%
      mutate(
        geno_prediction = as.integer(rowSums(
          across(all_of(selected_markers)), na.rm = TRUE
        ) > 0)
      )
  }

  # --- set up factors for yardstick ---
  # "1" as the first level = positive class
  df <- df %>%
    mutate(
      truth_value = factor(get(truth), levels = c(1, 0)),
      geno_prediction = factor(geno_prediction, levels = c(1, 0))
    )

  # --- compute yardstick metrics ---
  cm <- yardstick::conf_mat(df, truth = truth_value, estimate = geno_prediction)

  ys_metrics <- dplyr::bind_rows(
    yardstick::sens(df, truth = truth_value, estimate = geno_prediction),
    yardstick::spec(df, truth = truth_value, estimate = geno_prediction),
    yardstick::ppv(df, truth = truth_value, estimate = geno_prediction),
    yardstick::npv(df, truth = truth_value, estimate = geno_prediction),
    yardstick::accuracy(df, truth = truth_value, estimate = geno_prediction),
    yardstick::kap(df, truth = truth_value, estimate = geno_prediction),
    yardstick::f_meas(df, truth = truth_value, estimate = geno_prediction)
  )

  # --- compute AMR-specific metrics (VME and ME) ---
  sensitivity <- ys_metrics$.estimate[ys_metrics$.metric == "sens"]
  specificity <- ys_metrics$.estimate[ys_metrics$.metric == "spec"]

  vme <- 1 - sensitivity # FN / (TP + FN)
  me <- 1 - specificity  # FP / (TN + FP)

  amr_metrics <- tibble(
    .metric = c("VME", "ME"),
    .estimator = c("binary", "binary"),
    .estimate = c(vme, me)
  )

  all_metrics <- dplyr::bind_rows(ys_metrics, amr_metrics)

  # --- build return object ---
  # restore geno_prediction to integer in output data
  out_data <- binary_matrix %>%
    filter(!is.na(get(truth)))
  if (prediction_rule == "any") {
    out_data <- out_data %>%
      mutate(
        geno_prediction = as.integer(rowSums(
          across(all_of(selected_markers)), na.rm = TRUE
        ) > 0)
      )
  }

  result <- list(
    conf_mat = cm,
    metrics = all_metrics,
    data = out_data,
    markers_used = selected_markers,
    truth_col = truth,
    n = nrow(df)
  )

  structure(result, class = c("amr_concordance", class(result)))
}


#' Print method for amr_concordance objects
#'
#' Displays the confusion matrix and key concordance metrics.
#'
#' @param x An object of class `"amr_concordance"`.
#' @param ... Additional arguments (ignored).
#' @export
print.amr_concordance <- function(x, ...) {
  cat("AMR Genotype-Phenotype Concordance\n")
  cat(paste0("Truth: ", x$truth_col, " | Samples: ", x$n,
             " | Markers: ", length(x$markers_used), "\n"))
  cat(paste0("Markers used: ", paste(x$markers_used, collapse = ", "), "\n\n"))

  cat("Confusion Matrix:\n")
  print(x$conf_mat)
  cat("\n")

  # format and display key metrics
  m <- x$metrics
  fmt <- function(metric_name, digits = 4) {
    val <- m$.estimate[m$.metric == metric_name]
    if (length(val) == 0) return(NA_character_)
    format(round(val, digits), nsmall = digits)
  }

  cat("Metrics:\n")
  cat(paste0("  Sensitivity : ", fmt("sens"), "\n"))
  cat(paste0("  Specificity : ", fmt("spec"), "\n"))
  cat(paste0("  PPV         : ", fmt("ppv"), "\n"))
  cat(paste0("  NPV         : ", fmt("npv"), "\n"))
  cat(paste0("  Accuracy    : ", fmt("accuracy"), "\n"))
  cat(paste0("  Kappa       : ", fmt("kap"), "\n"))
  cat(paste0("  F-measure   : ", fmt("f_meas"), "\n"))
  cat(paste0("  VME         : ", fmt("VME"), "\n"))
  cat(paste0("  ME          : ", fmt("ME"), "\n"))

  invisible(x)
}
