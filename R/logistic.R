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

#' AMR Logistic Regression Analysis
#'
#' Performs logistic regression to analyze the relationship between genetic markers and phenotype (R, and NWT) for a specified antibiotic.
#' @param geno_table (Required if `binary_matrix` not provided) A data frame containing genotype data, formatted with [import_amrfp()]. Only used if `binary_matrix` not provided.
#' @param pheno_table (Required if `binary_matrix` not provided) A data frame containing phenotype data, formatted with [import_ast()]. Only used if `binary_matrix` not provided.
#' @param antibiotic (Required if `binary_matrix` not provided) A character string specifying the antibiotic of interest to filter phenotype data. The value must match one of the entries in the `drug_agent` column of `pheno_table`. Only used if `binary_matrix` not provided or if breakpoints required.
#' @param drug_class_list (Required if `binary_matrix` not provided) A character vector of drug classes to filter genotype data for markers related to the specified antibiotic. Markers in `geno_table` will be filtered based on whether their `drug_class` matches any value in this list. Only used if `binary_matrix` not provided.
#' @param geno_sample_col A character string (optional) specifying the column name in `geno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers. Only used if `binary_matrix` not provided.
#' @param pheno_sample_col A character string (optional) specifying the column name in `pheno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers. Only used if `binary_matrix` not provided.
#' @param sir_col A character string specifying the column name in `pheno_table` that contains the resistance interpretation (SIR) data. The values should be `"S"`, `"I"`, `"R"` or otherwise interpretable by [AMR::as.sir()]. If not provided, the first column prefixed with "phenotype*" will be used if present, otherwise an error is thrown.  Only used if `binary_matrix` not provided.
#' @param ecoff_col A character string specifying the column name in `pheno_table` that contains resistance interpretations (SIR) made against the ECOFF rather than a clinical breakpoint. The values should be `"S"`, `"I"`, `"R"` or otherwise interpretable by [AMR::as.sir()]. Default `ecoff`. Set to `NULL` if not available.  Only used if `binary_matrix` not provided.
#' @param marker_col (Optional) Name of the column containing the marker identifiers, whose unique values will be treated as predictors in the regression. Defaults to `"marker"`.
#' @param binary_matrix A data frame containing the original binary matrix output from the [get_binary_matrix()] function. If not provided (or set to `NULL`), user must specify `geno_table`, `pheno_table`, `antibiotic`, `drug_class_list` and optionally `geno_sample_col`, `pheno_sample_col`, `sir_col`, `ecoff_col`, `marker_col` to pass to [get_binary_matrix()].
#' @param maf (Optional) An integer specifying the minimum allele frequency (MAF) threshold. Markers with a MAF lower than this value will be excluded. Defaults to `10`.
#' @param fit_glm (Optional) Change to `TRUE` to fit model with glm. Otherwise fit model with logistf (default `FALSE`).
#' @param single_plot (Optional) A logical value. If `TRUE`, a single plot is produced comparing the estimates for resistance (`R`) and non-resistance (`NWT`). Otherwise, two plots are printed side-by-side. Defaults to `TRUE`.
#' @param colors (Optional) A vector of two colors, to use for R and NWT models in the plots. Defaults to `c("maroon", "blue4")`.
#' @param axis_label_size (Optional) A numeric value controlling the size of axis labels in the plot. Defaults to `9`.
#' @importFrom dplyr any_of select where
#' @importFrom ggplot2 ggtitle
#' @importFrom stats glm
#' @importFrom logistf logistf
#' @return A list with the following components:
#' - `binary_matrix`: The binary matrix of genetic data and phenotypic resistance information (either provided as input or generated by the function).
#' - `modelR`: The fitted logistic regression model for resistance (`R`).
#' - `modelNWT`: The fitted logistic regression model for non-resistance (`NWT`).
#' - `plot`: A ggplot object comparing the estimates for resistance and non-resistance with corresponding statistical significance indicators.
#' @export
#' @examples
#' # Example usage of the amr_logistic function
#' result <- amr_logistic(
#'   geno_table = import_amrfp(ecoli_geno_raw, "Name"),
#'   pheno_table = ecoli_ast,
#'   sir_col = "pheno_clsi",
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   maf = 10
#' )
#' # To access the plot:
#' print(result$plot)
amr_logistic <- function(geno_table, pheno_table,
                         antibiotic = NULL, drug_class_list = NULL,
                         geno_sample_col = NULL, pheno_sample_col = NULL,
                         sir_col = "pheno", ecoff_col = "ecoff",
                         marker_col = "marker.label", binary_matrix = NULL,
                         maf = 10, fit_glm = FALSE, single_plot = TRUE,
                         colors = c("maroon", "blue4"),
                         axis_label_size = 9) {
  # get binary matrix
  if (is.null(binary_matrix)) {
    cat("Generating geno-pheno binary matrix\n")
    binary_matrix <- get_binary_matrix(
      geno_table = geno_table,
      pheno_table = pheno_table,
      antibiotic = antibiotic,
      drug_class_list = drug_class_list,
      geno_sample_col = geno_sample_col,
      pheno_sample_col = pheno_sample_col,
      sir_col = sir_col,
      ecoff_col = ecoff_col,
      keep_assay_values = TRUE,
      marker_col = marker_col
    )
  }

  if (fit_glm) {
    cat("...Fitting logistic regression model to R using glm\n")
    if (sum(!is.na(binary_matrix$R)) > 0) {
      to_fit <- binary_matrix %>%
        select(-any_of(c("id", "pheno", "ecoff", "mic", "disk", "NWT"))) %>%
        select(R, where(~ sum(., na.rm = TRUE) >= maf))
      summarise_model_input(to_fit)
      modelR <- glm(R ~ ., data = to_fit, family = stats::binomial(link = "logit"))
      modelR <- glm_details(modelR) %>%
        mutate(marker = gsub("\\.\\.", ":", marker)) %>%
        mutate(marker = gsub("`", "", marker))
    }
    cat("...Fitting logistic regression model to NWT using glm\n")
    if (sum(!is.na(binary_matrix$NWT)) > 0) {
      to_fit <- binary_matrix %>%
        select(-any_of(c("id", "pheno", "ecoff", "mic", "disk", "R"))) %>%
        select(NWT, where(~ sum(., na.rm = TRUE) >= maf))
      summarise_model_input(to_fit)
      modelNWT <- glm(NWT ~ ., data = to_fit, family = stats::binomial(link = "logit"))
      modelNWT <- glm_details(modelNWT) %>%
        mutate(marker = gsub("\\.\\.", ":", marker)) %>%
        mutate(marker = gsub("`", "", marker))
    }
  } else {
    cat("...Fitting logistic regression model to R using logistf\n")
    if (sum(!is.na(binary_matrix$R)) > 0) {
      to_fit <- binary_matrix %>%
        filter(!is.na(R)) %>%
        select(-any_of(c("id", "pheno", "ecoff", "mic", "disk", "NWT"))) %>%
        select(R, where(~ sum(., na.rm = TRUE) >= maf))
      summarise_model_input(to_fit)
      modelR <- logistf::logistf(R ~ ., data = to_fit, pl = FALSE)
      modelR <- logistf_details(modelR) %>%
        mutate(marker = gsub("\\.\\.", ":", marker)) %>%
        mutate(marker = gsub("`", "", marker))
    }
    cat("...Fitting logistic regression model to NWT using logistf\n")
    if (sum(!is.na(binary_matrix$NWT)) > 0) {
      to_fit <- binary_matrix %>%
        filter(!is.na(NWT)) %>%
        select(-any_of(c("id", "pheno", "ecoff", "mic", "disk", "R"))) %>%
        select(NWT, where(~ sum(., na.rm = TRUE) >= maf))
      summarise_model_input(to_fit)
      modelNWT <- logistf::logistf(NWT ~ ., data = to_fit, pl = FALSE)
      modelNWT <- logistf_details(modelNWT) %>%
        mutate(marker = gsub("\\.\\.", ":", marker)) %>%
        mutate(marker = gsub("`", "", marker))
    }
  }

  cat("Generating plots\n")
  if (exists("modelR") & exists("modelNWT")) { # if we have 2 models, plot them together
    plot <- compare_estimates(modelR, modelNWT,
      single_plot = single_plot,
      title1 = "R", title2 = "NWT",
      colors = colors, axis_label_size = axis_label_size
    )
    if (single_plot) {
      label <- "Effect estimates for R and NWT"
      if (!is.null(antibiotic)) {
        label <- paste(label, "for", antibiotic)
      }
      if (!is.null(drug_class_list)) {
        subtitle <- paste("for", paste(drug_class_list, collapse = ","), "markers present in at least", maf, "samples")
      } else {
        subtitle <- paste("for markers present in at least", maf, "samples")
      }
      plot <- plot + ggtitle(label = label, subtitle = subtitle)
    }
  } else if (exists("modelR")) {
    plot <- plot_estimates(modelR)
  } else if (exists("modelNWT")) {
    plot <- plot_estimates(modelNWT)
  }
  if (!exists("modelNWT")) {
    modelNWT <- NULL
  } # need an object to return, set to null
  if (!exists("modelR")) {
    modelR <- NULL
  } # need an object to return, set to null

  print(plot)

  return(list(
    binary_matrix = binary_matrix,
    modelR = modelR,
    modelNWT = modelNWT,
    plot = plot
  ))
}


summarise_model_input <- function(dat) {
  cat(paste0(
    "   Filtered data contains ",
    nrow(dat),
    " samples (",
    sum(dat[, 1] == 1), " => 1, ",
    sum(dat[, 1] == 0), " => 0) and ",
    ncol(dat) - 1, " variables.\n"
  ))
}
