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
#' (R vs S, and/or NWT vs NWT) using a binary matrix from [get_binary_matrix()].
#' A genotypic prediction variable is defined on the basis of
#' genotype marker data (based on a variety of possible rules including
#' any/all/minimum number of markers; specifically those
#' markers or combinations exceeding a threshold positive
#' predictive value (PPV);  predictions from a logistic regression model; or a
#' user-defined field providing predictions).
#' This genotypic prediction is then
#' compared to the observed phenotypes using standard classification metrics
#' (via the `yardstick` pkg) and AMR-specific error rates (major error, ME
#' and very major error, VME) per ISO 20776-2 (and see
#' [FDA definitions](https://www.fda.gov/medical-devices/guidance-documents-medical-devices-and-radiation-emitting-products/antimicrobial-susceptibility-test-ast-systems-class-ii-special-controls-guidance-industry-and-fda)).
#' Supports evaluating both R and NWT outcomes
#' in a single call, with flexible prediction rules and marker inclusion options.
#'
#' @param binary_matrix A data frame output by [get_binary_matrix()], containing
#'   one row per sample, columns indicating binary phenotypes (`R`, `I`, `NWT`)
#'   and binary marker presence/absence.
#' @param prediction_rule The rule for generating genotypic predictions:
#'   `"any"` (default) predicts positive if any marker is present (after applying filters as specified by `markers`, `exclude_markers`, `ppv_threshold`, `pval_threshold`, `min_count`);
#'   `"all"` predicts positive only if all markers are present (after applying filters as specified by `markers`, `exclude_markers`, `ppv_threshold`, `pval_threshold`, `min_count`);
#'   a positive integer predicts positive if at least that many markers are present (after applying filters as specified by `markers`, `exclude_markers`, `ppv_threshold`, `pval_threshold`, `min_count`);
#'   `"logistic"` uses a logistic regression model from `logreg_results` to predict outcomes for each sample;
#'   `"combo_ppv"` predicts positive if a sample's marker combination has PPV >= `ppv_threshold` in `ppv_results`.
#' @param markers A character vector of marker column names to include in a
#'   new binary outcome variable `genotypic_prediction`.
#'   Default `NULL` includes all marker columns.
#' @param exclude_markers A character vector of marker column names to exclude
#'   from the genotypic prediction. Applied after `markers` filtering.
#' @param min_count An integer or `NULL`. Exclude markers with total frequency
#'   (column sum in binary_matrix) below this value. Default `NULL` (no filtering).
#' @param ppv_threshold A numeric PPV threshold (0-1). Used for solo PPV-based
#'   marker filtering when `solo_ppv_results` is provided, or as the combination
#'   PPV threshold when `prediction_rule = "combo_ppv"`.
#' @param solo_ppv_results Output of [solo_ppv_analysis()], used for PPV-based
#'   marker filtering when `ppv_threshold` is set.
#' @param ppv_results Output of [ppv()], required when
#'   `prediction_rule = "combo_ppv"`. The `summary` table from this object is
#'   used to identify marker combinations with PPV >= `ppv_threshold` for the
#'   relevant outcome (`R.ppv` or `NWT.ppv`). Samples whose marker combination
#'   matches any passing combination are predicted positive.
#' @param logreg_results Output of [amr_logistic()]. Used for p-value filtering
#'   (when `pval_threshold` is set) and for `prediction_rule = "logistic"`.
#' @param pval_threshold A numeric p-value threshold. Exclude markers with
#'   logistic regression p-value >= this value. Requires `logreg_results`.
#' @param prediction_col A character string naming a column in `binary_matrix`
#'   that contains a user-defined prediction (coded 0/1). When supplied, all
#'   marker filtering and prediction generation are bypassed; the specified column
#'   is used directly as the prediction for all outcomes in `truth`. This allows
#'   arbitrary prediction logic to be evaluated using the concordance metrics.
#' @param truth A character vector specifying the phenotypic truth column(s) to
#'   evaluate: `"R"` (resistant vs susceptible/intermediate), `"NWT"`
#'   (non-wildtype vs wildtype), or `c("R", "NWT")` (default) to evaluate both.
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
#' When `prediction_col` is supplied, all marker filtering and prediction
#' generation are bypassed. The named column (0/1-coded) is used directly as
#' the genotypic prediction. This is useful when the user has computed a
#' custom prediction (e.g. based on marker combinations from [ppv()]) and
#' wants to evaluate concordance metrics against the truth columns.
#'
#' When `prediction_rule = "combo_ppv"`, the function calls [get_combo_matrix()]
#' internally to derive a combination identifier for each sample. Samples whose
#' combination identifier matches any entry in `ppv_results$summary` with
#' outcome PPV >= `ppv_threshold` are predicted positive. The unique individual
#' markers contributing to passing combinations are reported in `markers_used`.
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
#' @importFrom dplyr across all_of any_of bind_rows filter mutate pull select
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
#' # Predict based on marker combinations with PPV >= 0.5 (from ppv())
#' ppv_res <- ppv(binary_matrix = binary_matrix)
#' result <- concordance(
#'   binary_matrix,
#'   prediction_rule = "combo_ppv",
#'   ppv_results = ppv_res,
#'   ppv_threshold = 0.5
#' )
#'
#' # Use a custom user-defined prediction column
#' binary_matrix$my_pred <- as.integer(binary_matrix$gyrA_S83L == 1 | binary_matrix$gyrA_D87N == 1)
#' result <- concordance(binary_matrix, prediction_col = "my_pred", truth = "R")
#'
#' # Access components
#' result$conf_mat$R
#' result$metrics
#' result$markers_used
#'
#' # predictions vs observed SIR calls
#' result$data %>% count(R_pred, pheno)
#' }
concordance <- function(binary_matrix,
                        markers = NULL,
                        exclude_markers = NULL,
                        ppv_threshold = NULL,
                        solo_ppv_results = NULL,
                        ppv_results = NULL,
                        prediction_col = NULL,
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
      stop("`prediction_rule` must be a positive integer, \"any\", \"all\", \"logistic\", or \"combo_ppv\".")
    }
  } else if (is.character(prediction_rule)) {
    prediction_rule <- match.arg(prediction_rule, choices = c("any", "all", "logistic", "combo_ppv"))
  } else {
    stop("`prediction_rule` must be a positive integer, \"any\", \"all\", \"logistic\", or \"combo_ppv\".")
  }

  if (identical(prediction_rule, "logistic") && is.null(logreg_results)) {
    stop("`logreg_results` must be provided when prediction_rule = \"logistic\".")
  }

  if (identical(prediction_rule, "combo_ppv") && is.null(ppv_results)) {
    stop("`ppv_results` must be provided when prediction_rule = \"combo_ppv\".")
  }

  if (identical(prediction_rule, "combo_ppv") && is.null(ppv_threshold)) {
    stop("`ppv_threshold` must be provided when prediction_rule = \"combo_ppv\".")
  }

  if (!is.null(prediction_col)) {
    if (!is.character(prediction_col) || length(prediction_col) != 1) {
      stop("`prediction_col` must be a single character string (column name).")
    }
    if (!(prediction_col %in% colnames(binary_matrix))) {
      stop(paste0("`prediction_col` column '", prediction_col, "' not found in binary_matrix."))
    }
  }

  if (!is.null(ppv_threshold) && is.null(solo_ppv_results) && !identical(prediction_rule, "combo_ppv")) {
    stop("`solo_ppv_results` must be provided when `ppv_threshold` is set (or use prediction_rule = \"combo_ppv\" with `ppv_results`).")
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

  if (is.null(prediction_col) && length(shared_markers) == 0) {
    stop("No markers remaining after filtering. Adjust `markers`, `exclude_markers`, or `min_count`.")
  }

  # --- pre-compute combo matrix for combo_ppv prediction rule ---
  if (identical(prediction_rule, "combo_ppv")) {
    combo_bm <- get_combo_matrix(binary_matrix)
  }

  # --- per-outcome loop ---
  conf_mat_list <- list()
  metrics_list <- list()
  markers_used_list <- list()
  n_vec <- integer(0)
  out_data <- binary_matrix

  for (outcome in truth) {
    # --- outcome-specific marker filtering (steps 4-5) ---
    # skipped when prediction_col is supplied or when using combo_ppv rule
    if (!is.null(prediction_col)) {
      # bypass all marker filtering; user-supplied column used directly
      selected_markers <- character(0)
      markers_used_list[[outcome]] <- character(0)
    } else if (identical(prediction_rule, "combo_ppv")) {
      # no per-marker filtering; combination PPV logic applied at prediction step
      selected_markers <- character(0)
    } else {
      selected_markers <- shared_markers

      # 4. PPV threshold filter (category-specific, solo markers)
      if (!is.null(ppv_threshold) && !is.null(solo_ppv_results)) {
        ppv_data <- solo_ppv_results$solo_stats
        passing_markers <- ppv_data %>%
          filter(category == outcome, ppv >= ppv_threshold) %>%
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
    }

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

    if (!is.null(prediction_col)) {
      # use user-supplied prediction column directly
      df[[pred_col]] <- df[[prediction_col]]
    } else if (identical(prediction_rule, "combo_ppv")) {
      # predict positive if sample's marker combination has outcome PPV >= threshold
      ppv_col_name <- paste0(outcome, ".ppv")
      if (!(ppv_col_name %in% colnames(ppv_results$summary))) {
        warning(paste0(
          "Column '", ppv_col_name, "' not found in ppv_results$summary for outcome '",
          outcome, "'. Skipping."
        ))
        next
      }
      passing_combos <- ppv_results$summary %>%
        filter(.data[[ppv_col_name]] >= ppv_threshold) %>%
        pull(combination_id)
      if (length(passing_combos) == 0) {
        warning(paste0(
          "No marker combinations pass ppv_threshold for outcome '", outcome, "'. Skipping."
        ))
        next
      }
      combo_df <- combo_bm %>% filter(!is.na(get(outcome)))
      df[[pred_col]] <- as.integer(combo_df$combination_id %in% passing_combos)
      # report unique individual markers from passing combinations
      markers_used_list[[outcome]] <- ppv_results$summary %>%
        filter(.data[[ppv_col_name]] >= ppv_threshold) %>%
        pull(marker_list) %>%
        strsplit(", ") %>%
        unlist() %>%
        unique()
    } else if (identical(prediction_rule, "logistic")) {
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
    if (!is.null(prediction_col)) {
      out_data[[pred_col]] <- binary_matrix[[prediction_col]]
    } else if (identical(prediction_rule, "combo_ppv")) {
      out_data[[pred_col]] <- as.integer(combo_bm$combination_id %in% passing_combos)
    } else if (identical(prediction_rule, "logistic")) {
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
  rule_label <- if (!is.null(prediction_col)) {
    paste0("user-supplied column: ", prediction_col)
  } else if (is.numeric(prediction_rule)) {
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
    mu <- x$markers_used[[outcome]]
    if (length(mu) == 0) {
      cat(paste0("Samples: ", x$n[outcome], "\n\n"))
    } else {
      cat(paste0(
        "Samples: ", x$n[outcome],
        " | Markers: ", length(mu), "\n"
      ))
      cat(paste0(
        "Markers used: ",
        paste(gsub("\\.\\.", ":", mu), collapse = ", "),
        "\n\n"
      ))
    }

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

#' Assess concordance between tables of observed and predicted phenotypes
#'
#' Calculates concordance from two data frames in standard phenotype table format, one with observed phenotype data and the other with predicted phenotype calls (e.g. output from AMRrules).
#'
#' @param pheno_table A tibble or data frame containing real (observed) phenotype data, in the format output by [import_pheno()].
#' @param pheno_pred_table A tibble or data frame containing predicted phenotype calls, e.g. AMRrules predictions imported in AMRgen using [import_amrrules_predictions()].
#' @param true_SIR_col Character. Name of the column containing S/I/R calls interpreted from real data against clinical breakpoints (default `"pheno_eucast"`).
#' @param true_ecoff_col Character. Name of the column containing WT/NWT calls interpreted from real data against ECOFF (default `"ecoff"`).
#' @param pred_SIR Character. Name of the column containing S/I/R calls predicted from genotypes (default `"clinical category"`).
#' @param pred_ecoff Character. Name of the column containing WT/NWT calls predicted from genotypes (default `"phenotype"`).
#' @param sample_col Character. Name of the column containing sample identifiers. This must be the same in both tables (default `"id"`).
#' @param drug_col Character. Name of the column containing drug agent identifiers. This must be the same in both tables (default `"drug_agent"`).
#' @param measure_col Character. Name of the column containing observed MIC or disk measurements for plotting. Valid options are `"mic"` or `"disk"` (default `"mic"`).
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{obs_pred}{A copy of the pheno_table with the predictions from pheno_pred_table merged in.}
#'   \item{metrics}{A tibble listing concordance metrics, comparing observed vs predicted R and/or NWT calls, for all drugs with matched samples in both input tables.}
#'   \item{drugs_R}{A long form tibble listing the number of samples with each combination of observed and predicted binary calls for R, for all drugs with matched samples in both input tables.}
#'   \item{drugs_NWT}{A long form tibble listing the number of samples with each combination of observed and predicted binary calls for NWT, for all drugs with matched samples in both input tables.}
#'   \item{plot_R}{A list containing two-panel plots summarising observed vs predicted clinical categories, one for each drug.}
#'   \item{plot_NWT}{A list containing two-panel plots summarising observed vs predicted WT/NWT categories, one for each drug.}
#'   \item{plot_R_dist}{A list containing plots of the assay value distributions (MIC or disk zones), coloured by S/I/R prediction, one for each drug.}
#'   \item{plot_NWT_dist}{A list containing plots of the assay value distributions (MIC or disk zones), coloured by WT/NWT prediction, one for each drug.}
#' }
#'
#' @importFrom dplyr summarise select count mutate left_join filter
#' @importFrom ggplot2 ggplot geom_tile scale_fill_gradient
#' @importFrom patchwork plot_layout plot_annotation
#' @examples
#' \dontrun{
#' concordance_from_tables(true_phenotypes, predicted_phenotypes)
#' }
#' @export
concordance_from_tables <- function(pheno_table,
                                    pheno_pred_table,
                                    true_SIR_col = "pheno_eucast",
                                    true_ecoff_col = "ecoff",
                                    pred_SIR = "clinical category",
                                    pred_ecoff = "phenotype",
                                    sample_col = "id",
                                    drug_col = "drug_agent",
                                    measure_col = "mic") {
  # join observed & predicted pheno tables
  obs_pred <- pheno_table %>%
    left_join(pheno_pred_table, by = c(sample_col, drug_col))

  # encode binary outcomes from each
  obs_pred <- obs_pred %>%
    mutate(!!sym(pred_SIR) := if_else(is.na(!!sym(pred_SIR)), "S", !!sym(pred_SIR))) %>%
    mutate(!!sym(pred_ecoff) := if_else(is.na(!!sym(pred_ecoff)), "NWT", !!sym(pred_ecoff))) %>%
    mutate(R = case_when(
      !!sym(true_SIR_col) == "R" ~ 1,
      !!sym(true_SIR_col) == "I" ~ 0,
      !!sym(true_SIR_col) == "S" ~ 0,
      TRUE ~ NA
    )) %>%
    mutate(NWT = case_when(
      !!sym(true_ecoff_col) == "NWT" ~ 1,
      !!sym(true_ecoff_col) == "WT" ~ 0,
      TRUE ~ NA
    )) %>%
    mutate(predR = case_when(
      !!sym(pred_SIR) == "R" ~ 1,
      !!sym(pred_SIR) == "S" ~ 0,
      is.na(!!sym(pred_SIR)) ~ 0,
      TRUE ~ NA
    )) %>%
    mutate(predNWT = case_when(
      !!sym(pred_ecoff) == "nonwildtype" ~ 1,
      !!sym(pred_ecoff) == "wildtype" ~ 0,
      is.na(!!sym(pred_ecoff)) ~ 0,
      TRUE ~ NA
    ))

  drugs_R <- obs_pred %>%
    filter(!is.na(R) & !is.na(predR)) %>%
    count(!!sym(drug_col), R, predR)

  drugs_NWT <- obs_pred %>%
    filter(!is.na(NWT) & !is.na(predNWT)) %>%
    count(!!sym(drug_col), NWT, predNWT)

  plot_R <- list()
  plot_NWT <- list()
  dist_plot_R <- list()
  dist_plot_NWT <- list()
  concordance_metrics <- tibble()

  drug_list <- unique(c(drugs_R[[drug_col]], drugs_NWT[[drug_col]]))

  for (this_drug in as.character(drug_list)) {
    drug_data <- obs_pred %>%
      filter(!!sym(drug_col) == this_drug)

    # plot distribution of observed S/I/R values for genomes called as S (incl NA) or R
    bar_plot_R <- NULL
    bar_plot_R <- drug_data %>%
      filter(!is.na(!!sym(true_SIR_col)) & !is.na(!!sym(true_SIR_col))) %>%
      count(!!sym(true_SIR_col), !!sym(pred_SIR)) %>%
      ggplot(aes(y = !!sym(pred_SIR), x = n, fill = !!sym(true_SIR_col))) +
      geom_col() +
      theme_bw() +
      labs(y = "Predicted", fill = "Observed", x = "Count")

    crosstabs_plot_R <- NULL
    crosstabs_plot_R <- drug_data %>%
      count(!!sym(true_SIR_col), !!sym(pred_SIR)) %>%
      group_by(!!sym(pred_SIR)) %>%
      mutate(row_pct = n / sum(n)) %>%
      ggplot(aes(x = !!sym(true_SIR_col), y = !!sym(pred_SIR), fill = n)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%d\n%.1f%%", n, 100 * row_pct))) +
      scale_fill_gradient(low = "white", high = "skyblue") +
      labs(y = "Predicted", fill = "Counts", x = "Observed") +
      theme_bw()

    plot_R[[this_drug]] <- crosstabs_plot_R / bar_plot_R +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = ab_name(as.ab(this_drug)),
        subtitle = "observed vs predicted clinical phenotype"
      )

    if (measure_col %in% colnames(obs_pred)) {
      dist_plot_R[[this_drug]] <- obs_pred %>%
        assay_by_var(
          antibiotic = this_drug, colour_by = pred_SIR, measure = measure_col,
          plot_title = paste(
            ab_name(as.ab(this_drug)), measure_col,
            "distribution vs predicted clinical phenotype"
          )
        )
    }

    # plot distribution of observed S/I/R values for genomes called as S (incl NA) or R
    bar_plot_NWT <- NULL
    bar_plot_NWT <- drug_data %>%
      filter(!is.na(!!sym(true_ecoff_col)) & !is.na(!!sym(pred_ecoff))) %>%
      count(!!sym(true_ecoff_col), !!sym(pred_ecoff)) %>%
      ggplot(aes(y = !!sym(pred_ecoff), x = n, fill = !!sym(true_ecoff_col))) +
      geom_col() +
      labs(y = "Predicted", fill = "Observed", x = "Count")

    crosstabs_plot_NWT <- NULL
    crosstabs_plot_NWT <- drug_data %>%
      count(!!sym(true_ecoff_col), !!sym(pred_ecoff)) %>%
      group_by(!!sym(pred_ecoff)) %>%
      mutate(row_pct = n / sum(n)) %>%
      ggplot(aes(x = !!sym(true_ecoff_col), y = !!sym(pred_ecoff), fill = n)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%d\n%.1f%%", n, 100 * row_pct))) +
      scale_fill_gradient(low = "white", high = "skyblue") +
      labs(y = "Predicted", fill = "Counts", x = "Observed") +
      theme_bw()

    plot_NWT[[this_drug]] <- crosstabs_plot_NWT / bar_plot_NWT +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = ab_name(as.ab(this_drug)),
        subtitle = "observed vs predicted WT/NWT phenotype"
      )

    if (measure_col %in% colnames(obs_pred)) {
      dist_plot_NWT[[this_drug]] <- obs_pred %>%
        assay_by_var(
          antibiotic = this_drug, colour_by = pred_ecoff, measure = measure_col,
          plot_title = paste(
            ab_name(as.ab(this_drug)),
            measure_col, "distribution vs predicted WT/NWT phenotype"
          )
        )
    }

    conR <- NULL
    if (nrow(drugs_R) > 0) {
      conR <- safe_execute(concordance(drug_data %>% filter(!is.na(R) & !is.na(predR)),
        truth = "R", prediction_col = "predR"
      )$metrics)
      if (!is.null(conR)) {
        concordance_metrics <- bind_rows(concordance_metrics, conR %>% mutate(drug_agent = this_drug)) %>%
          relocate(drug_agent)
      }
    }

    conNWT <- NULL
    if (nrow(drugs_NWT) > 0) {
      conNWT <- safe_execute(concordance(drug_data %>% filter(!is.na(NWT) & !is.na(predNWT)),
        truth = "NWT", prediction_col = "predNWT"
      )$metrics)
      if (!is.null(conNWT)) {
        concordance_metrics <- bind_rows(concordance_metrics, conNWT %>% mutate(drug_agent = this_drug)) %>%
          relocate(drug_agent)
      }
    }
  }

  return(list(
    obs_pred = obs_pred,
    metrics = concordance_metrics,
    drugs_R = drugs_R,
    drugs_NWT = drugs_NWT,
    plot_R = plot_R,
    plot_NWT = plot_NWT,
    plot_R_dist = dist_plot_R,
    plot_NWT_dist = dist_plot_NWT
  ))
}
