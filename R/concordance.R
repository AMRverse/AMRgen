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

#' Calculate genotype-phenotype concordance from binary matrix
#'
#' Compares genotypes (presence of resistance markers) to observed phenotypes
#' (resistant vs susceptible) using a binary matrix from [get_binary_matrix()].
#' A genotypic prediction variable is defined on the basis of presence of
#' genotype markers (either any marker in the input table, or those defined by
#' an input inclusion list or exclusion list). This genotypic prediction is then
#' compared to the observed phenotypes using standard classification metrics
#' (via the `yardstick` pkg) and AMR-specific error rates (major error, ME
#' and very major error, VME) per ISO 20776-2 (and see
#' [FDA definitions](https://www.fda.gov/medical-devices/guidance-documents-medical-devices-and-radiation-emitting-products/antimicrobial-susceptibility-test-ast-systems-class-ii-special-controls-guidance-industry-and-fda).
#' Supports evaluating both R and NWT outcomes
#' in a single call, with flexible prediction rules and marker inclusion options.
#'
#' @param binary_matrix A data frame output by [get_binary_matrix()], containing
#'   one row per sample, columns indicating binary phenotypes (`R`, `I`, `NWT`)
#'   and binary marker presence/absence.
#' @param markers A character vector of marker column names to include in a
#'   new binary outcome variable `genotypic_prediction`.
#'   Default `NULL` includes all marker columns.
#' @param exclude_markers A character vector of marker column names to exclude
#'   from the genotypic prediction. Applied after `markers` filtering.
#' @param ppv_threshold A numeric PPV threshold (0-1). Markers with solo PPV
#'   below this value are excluded. Requires `solo_ppv_results`.
#' @param solo_ppv_results Output of [solo_ppv_analysis()], used for PPV-based
#'   marker filtering when `ppv_threshold` is set.
#' @param truth A character vector specifying the phenotypic truth column(s) to
#'   evaluate: `"R"` (resistant vs susceptible/intermediate), `"NWT"`
#'   (non-wildtype vs wildtype), or `c("R", "NWT")` (default) to evaluate both.
#' @param prediction_rule The rule for generating genotypic predictions:
#'   `"any"` (default) predicts positive if any marker is present (after applying filters as specified by `markers`, `exclude_markers`, `ppv_threshold`, `pval_threshold`, `min_count`);
#'   `"all"` predicts positive only if all markers are present (after applying filters as specified by `markers`, `exclude_markers`, `ppv_threshold`, `pval_threshold`, `min_count`);
#'   a positive integer predicts positive if at least that many markers are present (after applying filters as specified by `markers`, `exclude_markers`, `ppv_threshold`, `pval_threshold`, `min_count`);
#'   `"logistic"` uses a logistic regression model from `logreg_results` to predict outcomes for each sample.
#' @param min_count An integer or `NULL`. Exclude markers with total frequency
#'   (column sum in binary_matrix) below this value. Default `NULL` (no filtering).
#' @param logreg_results Output of [amr_logistic()]. Used for p-value filtering
#'   (when `pval_threshold` is set) and for `prediction_rule = "logistic"`.
#' @param pval_threshold A numeric p-value threshold. Exclude markers with
#'   logistic regression p-value >= this value. Requires `logreg_results`.
#'
#' @details
#' The function identifies marker columns as all columns not in the reserved set
#' (`id`, `pheno`, `ecoff`, `R`, `I`, `NWT`, `mic`, `disk`). It then applies
#' filtering in order: inclusion list `markers`, exclusion list `exclude_markers`, `min_count`,
#' `ppv_threshold` filtering, then `pval_threshold` filtering.
#'
#' Marker filtering is performed per outcome (R, NWT), since PPV-based and
#' p-value-based filters are category/outcome-specific.
#'
#' Standard metrics (sensitivity, specificity, PPV, NPV, accuracy, kappa,
#' F-measure) are calculated using pkg `yardstick`. AMR-specific error rates are
#' computed internally:
#' - **VME** (Very Major Error): FN / (TP + FN) = 1 - sensitivity. Proportion of
#'   truly resistant isolates not predicted as such from genotype.
#' - **ME** (Major Error): FP / (TN + FP) = 1 - specificity. Proportion of
#'   truly susceptible isolates incorrectly predicted resistant from genotype.
#'
#' @return An S3 object of class `"amr_concordance"`, a list containing:
#' - `conf_mat`: Named list of yardstick confusion matrix objects (e.g.
#'   `list(R = <cm>, NWT = <cm>)`).
#' - `metrics`: A tibble with columns `outcome`, `metric`, `estimate`.
#' - `data`: The input binary matrix with added `R_pred` and/or `NWT_pred` columns.
#' - `markers_used`: Named list of character vectors of markers used per outcome.
#' - `truth_col`: Character vector of truth columns evaluated.
#' - `n`: Named integer vector of sample counts per outcome.
#' - `prediction_rule`: The prediction rule used.
#'
#' @importFrom dplyr across all_of any_of bind_rows filter mutate select
#' @importFrom stats predict
#' @importFrom tibble tibble
#' @importFrom yardstick conf_mat sens spec ppv npv accuracy kap f_meas
#' @seealso [get_binary_matrix()], [solo_ppv_analysis()], [amr_logistic()], [yardstick]
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
#' # Basic concordance for both R and NWT
#' result <- concordance(binary_matrix)
#' result
#'
#' # Single outcome
#' result <- concordance(binary_matrix, truth = "R")
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
#' # Require at least 2 markers present for prediction
#' result <- concordance(binary_matrix, prediction_rule = 2)
#'
#' # Use logistic regression model for prediction
#' logreg <- amr_logistic(binary_matrix = binary_matrix)
#' result <- concordance(
#'   binary_matrix,
#'   prediction_rule = "logistic",
#'   logreg_results = logreg
#' )
#'
#' # Access components
#' result$conf_mat$R
#' result$metrics
#' result$markers_used
#' }
concordance <- function(binary_matrix,
                        markers = NULL,
                        exclude_markers = NULL,
                        ppv_threshold = NULL,
                        solo_ppv_results = NULL,
                        truth = c("R", "NWT"),
                        prediction_rule = "any",
                        min_count = NULL,
                        logreg_results = NULL,
                        pval_threshold = NULL) {
  # --- input validation ---
  if (!is.data.frame(binary_matrix)) {
    stop("`binary_matrix` must be a data frame (output of get_binary_matrix()).")
  }

  truth <- match.arg(truth, choices = c("R", "NWT"), several.ok = TRUE)

  for (tc in truth) {
    if (!(tc %in% colnames(binary_matrix))) {
      stop(paste0("Truth column '", tc, "' not found in binary_matrix."))
    }
  }

  # validate prediction_rule
  if (is.numeric(prediction_rule)) {
    if (prediction_rule < 1 || prediction_rule != as.integer(prediction_rule)) {
      stop("`prediction_rule` must be a positive integer, \"any\", \"all\", or \"logistic\".")
    }
  } else if (is.character(prediction_rule)) {
    prediction_rule <- match.arg(prediction_rule, choices = c("any", "all", "logistic"))
  } else {
    stop("`prediction_rule` must be a positive integer, \"any\", \"all\", or \"logistic\".")
  }

  if (identical(prediction_rule, "logistic") && is.null(logreg_results)) {
    stop("`logreg_results` must be provided when prediction_rule = \"logistic\".")
  }

  if (!is.null(ppv_threshold) && is.null(solo_ppv_results)) {
    stop("`solo_ppv_results` must be provided when `ppv_threshold` is set.")
  }

  if (!is.null(pval_threshold) && is.null(logreg_results)) {
    stop("`logreg_results` must be provided when `pval_threshold` is set.")
  }

  if (!is.null(min_count) && (!is.numeric(min_count) || min_count < 1)) {
    stop("`min_count` must be a positive integer or NULL.")
  }

  # --- identify marker columns ---
  reserved_cols <- c("id", "pheno", "ecoff", "R", "I", "NWT", "mic", "disk")
  all_markers <- setdiff(colnames(binary_matrix), reserved_cols)

  if (length(all_markers) == 0) {
    stop("No marker columns found in binary_matrix.")
  }

  # --- shared marker filtering (steps 1-3, outcome-independent) ---
  shared_markers <- all_markers

  # 1. custom markers list
  if (!is.null(markers)) {
    unknown <- setdiff(markers, all_markers)
    if (length(unknown) > 0) {
      warning(paste(
        "Markers not found in binary_matrix (ignored):",
        paste(unknown, collapse = ", ")
      ))
    }
    shared_markers <- intersect(markers, all_markers)
  }

  # 2. exclude markers
  if (!is.null(exclude_markers)) {
    shared_markers <- setdiff(shared_markers, exclude_markers)
  }

  # 3. min_count filter
  if (!is.null(min_count)) {
    marker_sums <- colSums(binary_matrix[, shared_markers, drop = FALSE], na.rm = TRUE)
    shared_markers <- names(marker_sums[marker_sums >= min_count])
  }

  if (length(shared_markers) == 0) {
    stop("No markers remaining after filtering. Adjust `markers`, `exclude_markers`, or `min_count`.")
  }

  # --- per-outcome loop ---
  conf_mat_list <- list()
  metrics_list <- list()
  markers_used_list <- list()
  n_vec <- integer(0)
  out_data <- binary_matrix

  for (outcome in truth) {
    # outcome-specific marker filtering (steps 4-5)
    selected_markers <- shared_markers

    # 4. PPV threshold filter (category-specific)
    if (!is.null(ppv_threshold)) {
      ppv_category <- outcome
      ppv_data <- solo_ppv_results$solo_stats
      passing_markers <- ppv_data %>%
        filter(category == ppv_category, ppv >= ppv_threshold) %>%
        pull(marker)
      # normalize : to .. for matching binary_matrix column names
      passing_markers_norm <- gsub(":", "..", passing_markers)
      selected_markers <- intersect(selected_markers, passing_markers_norm)
    }

    # 5. pval threshold filter (outcome-specific)
    if (!is.null(pval_threshold)) {
      logreg_model_key <- paste0("model", outcome)
      logreg_summary <- logreg_results[[logreg_model_key]]
      if (!is.null(logreg_summary)) {
        passing_lr <- logreg_summary %>%
          filter(marker != "(Intercept)", pval < pval_threshold) %>%
          pull(marker)
        # normalize : to .. for matching binary_matrix column names
        passing_lr_norm <- gsub(":", "..", passing_lr)
        selected_markers <- intersect(selected_markers, passing_lr_norm)
      }
    }

    if (length(selected_markers) == 0) {
      warning(paste0(
        "No markers remaining for outcome '", outcome,
        "' after filtering. Skipping."
      ))
      next
    }

    markers_used_list[[outcome]] <- selected_markers

    # --- filter samples with non-NA truth ---
    df <- binary_matrix %>%
      filter(!is.na(get(outcome)))

    if (nrow(df) == 0) {
      warning(paste0(
        "No samples with non-NA values in truth column '",
        outcome, "'. Skipping."
      ))
      next
    }

    n_vec[outcome] <- nrow(df)

    # --- generate prediction ---
    pred_col <- paste0(outcome, "_pred")

    if (identical(prediction_rule, "logistic")) {
      # use raw logistic regression model to predict
      raw_model_key <- paste0("raw_model", outcome)
      raw_model <- logreg_results[[raw_model_key]]
      if (is.null(raw_model)) {
        stop(paste0(
          "No raw model found in logreg_results$", raw_model_key,
          ". Ensure amr_logistic() returned raw model objects."
        ))
      }
      probs <- predict(raw_model, newdata = df, type = "response")
      df[[pred_col]] <- as.integer(probs > 0.5)
    } else if (identical(prediction_rule, "any")) {
      df[[pred_col]] <- as.integer(
        rowSums(df[, selected_markers, drop = FALSE], na.rm = TRUE) > 0
      )
    } else if (identical(prediction_rule, "all")) {
      n_markers <- length(selected_markers)
      df[[pred_col]] <- as.integer(
        rowSums(df[, selected_markers, drop = FALSE], na.rm = TRUE) == n_markers
      )
    } else if (is.numeric(prediction_rule)) {
      df[[pred_col]] <- as.integer(
        rowSums(df[, selected_markers, drop = FALSE], na.rm = TRUE) >= prediction_rule
      )
    }

    # --- set up factors for yardstick ---
    df <- df %>%
      mutate(
        truth_value = factor(get(outcome), levels = c(1, 0)),
        geno_prediction = factor(get(pred_col), levels = c(1, 0))
      )

    # --- compute confusion matrix ---
    cm <- yardstick::conf_mat(df, truth = truth_value, estimate = geno_prediction)
    conf_mat_list[[outcome]] <- cm

    # --- compute yardstick metrics ---
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
    me <- 1 - specificity # FP / (TN + FP)

    amr_metrics <- tibble(
      .metric = c("VME", "ME"),
      .estimator = c("binary", "binary"),
      .estimate = c(vme, me)
    )

    # --- clean metrics format ---
    outcome_metrics <- dplyr::bind_rows(ys_metrics, amr_metrics) %>%
      mutate(outcome = outcome) %>%
      select(outcome, metric = .metric, estimate = .estimate)

    metrics_list[[outcome]] <- outcome_metrics

    # --- add prediction column to output data ---
    # compute prediction for the full binary_matrix (including NA truth rows)
    if (identical(prediction_rule, "logistic")) {
      raw_model_key <- paste0("raw_model", outcome)
      raw_model <- logreg_results[[raw_model_key]]
      probs_all <- predict(raw_model, newdata = binary_matrix, type = "response")
      out_data[[pred_col]] <- as.integer(probs_all > 0.5)
    } else if (identical(prediction_rule, "any")) {
      out_data[[pred_col]] <- as.integer(
        rowSums(binary_matrix[, selected_markers, drop = FALSE], na.rm = TRUE) > 0
      )
    } else if (identical(prediction_rule, "all")) {
      n_markers <- length(selected_markers)
      out_data[[pred_col]] <- as.integer(
        rowSums(binary_matrix[, selected_markers, drop = FALSE], na.rm = TRUE) == n_markers
      )
    } else if (is.numeric(prediction_rule)) {
      out_data[[pred_col]] <- as.integer(
        rowSums(binary_matrix[, selected_markers, drop = FALSE], na.rm = TRUE) >= prediction_rule
      )
    }
    out_data <- out_data %>% relocate(any_of(c("R", "NWT")), .after = 1)
    out_data <- out_data %>% relocate(any_of(c("R_pred", "NWT_pred")), .after = 1)
  }

  if (length(conf_mat_list) == 0) {
    stop("No outcomes could be evaluated. Check truth columns and filtering parameters.")
  }

  # --- combine metrics ---
  all_metrics <- dplyr::bind_rows(metrics_list)

  # --- format prediction_rule for display ---
  rule_label <- if (is.numeric(prediction_rule)) {
    as.character(prediction_rule)
  } else {
    prediction_rule
  }

  # --- build return object ---
  result <- list(
    conf_mat = conf_mat_list,
    metrics = all_metrics,
    data = out_data,
    markers_used = markers_used_list,
    truth_col = truth,
    n = n_vec,
    prediction_rule = rule_label
  )

  structure(result, class = c("amr_concordance", class(result)))
}


#' Print method for amr_concordance objects
#'
#' Displays the confusion matrix and key concordance metrics for each outcome.
#'
#' @param x An object of class `"amr_concordance"`.
#' @param ... Additional arguments (ignored).
#' @export
print.amr_concordance <- function(x, ...) {
  cat("AMR Genotype-Phenotype Concordance\n")
  cat(paste0("Prediction rule: ", x$prediction_rule, "\n"))

  for (outcome in names(x$conf_mat)) {
    cat(paste0("\n--- Outcome: ", outcome, " ---\n"))
    cat(paste0(
      "Samples: ", x$n[outcome],
      " | Markers: ", length(x$markers_used[[outcome]]), "\n"
    ))
    cat(paste0(
      "Markers used: ",
      paste(gsub("\\.\\.", ":", x$markers_used[[outcome]]), collapse = ", "),
      "\n\n"
    ))

    cat("Confusion Matrix:\n")
    print(x$conf_mat[[outcome]])
    cat("\n")

    # format and display metrics for this outcome
    m <- x$metrics[x$metrics$outcome == outcome, ]
    fmt <- function(metric_name, digits = 4) {
      val <- m$estimate[m$metric == metric_name]
      if (length(val) == 0) {
        return(NA_character_)
      }
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
  }

  invisible(x)
}
