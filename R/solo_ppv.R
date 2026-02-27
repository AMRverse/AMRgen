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

#' Perform Solo PPV Analysis for AMR Markers
#'
#' This function performs a Positive Predictive Value (PPV) analysis for AMR markers associated with a given antibiotic and drug class. It calculates the PPV for solo markers and visualizes the results using various plots.
#' @param binary_matrix A data frame containing the original binary matrix output from the [get_binary_matrix()] function. If not provided (or set to `NULL`), user must specify `geno_table`, `pheno_table`, `antibiotic`, `drug_class_list` and optionally `geno_sample_col`, `pheno_sample_col`, `sir_col`, `ecoff_col`, `marker_col` to pass to [get_binary_matrix()].
#' @param geno_table (Required if `binary_matrix` not provided) A data frame containing genotype data, formatted with [import_amrfp()]. Only used if `binary_matrix` not provided.
#' @param pheno_table (Required if `binary_matrix` not provided) A data frame containing phenotype data, formatted with [import_ast()]. Only used if `binary_matrix` not provided.
#' @param antibiotic (Required if `binary_matrix` not provided) A character string specifying the antibiotic of interest to filter phenotype data. The value must match one of the entries in the `drug_agent` column of `pheno_table`. Only used if `binary_matrix` not provided or if breakpoints required.
#' @param drug_class_list (Required if `binary_matrix` not provided) A character vector of drug classes to filter genotype data for markers related to the specified antibiotic. Markers in `geno_table` will be filtered based on whether their `drug_class` matches any value in this list. Only used if `binary_matrix` not provided.
#' @param geno_sample_col A character string (optional) specifying the column name in `geno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers. Only used if `binary_matrix` not provided.
#' @param pheno_sample_col A character string (optional) specifying the column name in `pheno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers. Only used if `binary_matrix` not provided.
#' @param sir_col A character string specifying the column name in `pheno_table` that contains the resistance interpretation (SIR) data. The values should be `"S"`, `"I"`, `"R"` or otherwise interpretable by [AMR::as.sir()]. If not provided, the first column prefixed with "phenotype*" will be used if present, otherwise an error is thrown.  Only used if `binary_matrix` not provided.
#' @param ecoff_col A character string specifying the column name in `pheno_table` that contains the ECOFF interpretation of phenotype. The values should be interpretable as `"WT"` (wildtype) and `"NWT"` (nonwildtype), or `"S"` / `"I"` / `"R"`. Default `ecoff`. Set to `NULL` if not available.  Only used if `binary_matrix` not provided.
#' @param marker_col A character string specifying the column name in `geno_table` containing the marker identifiers. Default `"marker"`. Only used if `binary_matrix` not provided.
#' @param reverse_order A logical indicating whether to reverse the order of rows in the plot, so that markers are ordered from lowest R PPV to highest (default `FALSE`, i.e. markers are ordered from highest to lowest PPV).
#' @param icat A logical indicating whether to calculate PPV for `"I"` (if such a category exists in the phenotype column) (default `FALSE`).
#' @param min Minimum number of genomes with the solo marker, to include the marker in the plot (default `1`).
#' @param pd A `ggplot2::position_dodge()` object controlling horizontal spacing of points and confidence intervals in the PPV plot. Default `position_dodge(width = 0.8)`.
#' @param axis_label_size Font size for axis labels in the PPV plot (default `9`).
#' @param excludeRanges Vector of phenotype categories (comprised of `"R"`, `"I"`, `"NWT"`) for which we should ignore MIC values expressed as ranges when calculating PPVs. To include MICs expressed as ranges set this to `NULL`.
#' @param colours_SIR A named vector of colours for the percentage bar plot. The names should be the phenotype categories (e.g., `"R"`, `"I"`, `"S"`), and the values should be valid colour names or hexadecimal colour codes. Default values are those used in the AMR package [scale_fill_sir()].
#' @param colours_ppv A named vector of colours for the plot of PPV estimates. The names should be `"R"`, `"I"`, `"NWT"`, and the values should be valid colour names or hexadecimal colour codes.
#' @details The function analyzes the predictive power of individual AMR markers when they are found 'solo' in the genome with no other markers associated with the same class. The phenotype data are matched with genotype presence/absence and then stratified to compute PPV for resistance and non-wild-type interpretations. The function also generates plots to aid in interpretation.
#' @importFrom AMR scale_fill_sir
#' @importFrom dplyr any_of bind_rows filter group_by mutate n relocate rename select summarise
#' @importFrom ggplot2 aes after_stat element_text geom_bar geom_linerange geom_point geom_text geom_vline ggplot ggtitle labs position_dodge position_fill scale_colour_manual scale_y_discrete theme theme_bw theme_light xlim
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom patchwork plot_layout
#' @importFrom tibble tibble
#' @return A list containing the following elements:
#' - `solo_stats`: A data frame summarizing the PPV for resistance (R vs S/I) and NWT (R/I vs S), including the number of positive hits, sample size, PPV, and 95% confidence intervals for each marker.
#' - `combined_plot`: A combined ggplot object showing the PPV plot for the solo markers, and a bar plot for the phenotype distribution.
#' - `solo_binary`: A dataframe with binary values indicating the presence or absence of the solo markers.
#' - `solo_binary_norange`: A dataframe with binary values indicating the presence or absence of the solo markers, excluding samples for which MIC or disk measures are expressed as ranges.
#' - `amr_binary`: A copy of the genotype-phenotype binary matrix for all markers (either provided as input or generated by the function)
#' - `plot_order`: Ordered list of rows in the output plot (to facilitate alignment with other plots)
#' @export
#' @examples
#' \dontrun{
#' geno_table <- import_amrfp(ecoli_geno_raw, "Name")
#' head(geno_table)
#'
#' # Generate binary matrix
#' binary_matrix <- get_binary_matrix(
#'   geno_table = geno_table,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno_clsi",
#'   keep_assay_values = TRUE,
#'   keep_assay_values_from = "mic"
#' )
#'
#' # Run solo PPV analysis plot analysis using this binary_matrix
#' # (note antibiotic and drug_class list are optional here, and only used
#' # for titling the plot)
#' soloPPV_cipro <- solo_ppv_analysis(binary_matrix = binary_matrix)
#'
#' soloPPV_cipro$solo_stats
#' soloPPV_cipro$combined_plot
#' }
solo_ppv_analysis <- function(geno_table, pheno_table,
                              antibiotic = NULL, drug_class_list = NULL,
                              geno_sample_col = NULL, pheno_sample_col = NULL,
                              sir_col = NULL, ecoff_col = "ecoff", icat = FALSE,
                              marker_col = "marker", reverse_order = FALSE,
                              binary_matrix = NULL,
                              min = 1,
                              axis_label_size = 9,
                              pd = position_dodge(width = 0.8),
                              excludeRanges = NULL,
                              colours_SIR = c(
                                S = "#3CAEA3", SDD = "#8FD6C4",
                                I = "#F6D55C", R = "#ED553B"
                              ),
                              colours_ppv = c(
                                "R" = "maroon", "I" = "skyblue",
                                "NWT" = "navy"
                              )) {
  # get binary matrix
  if (is.null(binary_matrix)) {
    cat("Generating geno-pheno binary matrix\n")

    # check there is a SIR column specified
    if (is.null(sir_col)) {
      # make a sensible guess
      sir_col <- pheno_table %>%
        select(starts_with("pheno")) %>%
        colnames() %>%
        first()
      if (!is.na(sir_col)) {
        cat(paste("WARNING: `sir_col` not provided, using first column with prefix 'pheno':", sir_col))
      } else {
        stop("`sir_col` not provided. Please specify a column with S/I/R phenotype values.")
      }
    }
    if (!(sir_col %in% colnames(pheno_table))) {
      stop(paste0("Column: '", sir_col, "' not found in input phenotype data. Please specify a valid column with S/I/R values."))
    }

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

  if (icat & ("I" %in% binary_matrix$pheno)) { ### TO DO: CHECK HOW TO DO THIS WITH COLNAME
    binary_matrix <- binary_matrix %>%
      mutate(I = case_when(
        pheno %in% c("I", "R") ~ 1,
        pheno == "S" ~ 0,
        TRUE ~ NA
      )) %>%
      relocate(I, .after = R)
  } else {
    icat <- FALSE
  }

  # get solo markers
  marker_counts <- binary_matrix %>%
    select(-any_of(c("id", "pheno", "ecoff", "R", "I", "NWT", "mic", "disk"))) %>%
    rowSums()

  solo_binary <- binary_matrix %>%
    filter(marker_counts == 1) %>%
    pivot_longer(!any_of(c("id", "pheno", "ecoff", "R", "I", "NWT", "solo", "mic", "disk")), names_to = "marker") %>%
    mutate(marker = gsub("\\.\\.", ":", marker)) %>%
    mutate(marker = gsub("`", "", marker)) %>%
    filter(value == 1) %>%
    filter(dplyr::if_any(any_of(c("ecoff", "pheno")), ~ !is.na(.)))

  if (nrow(solo_binary) == 0) {
    stop("No solo markers found")
  }

  # summarise numerator, denominator, proportion, 95% CI - for R and NWT

  if (sum(c("R", "I", "NWT") %in% excludeRanges) > 0 & "mic" %in% colnames(solo_binary)) {
    solo_binary_norange <- solo_binary %>% filter(!grepl("<", mic))
  } else {
    solo_binary_norange <- NULL
  }

  if (sum(!is.na(solo_binary$pheno)) > 0) {
    if ("R" %in% excludeRanges & "mic" %in% colnames(solo_binary)) {
      solo_binary_R <- solo_binary_norange
    } else {
      solo_binary_R <- solo_binary
    }
    solo_stats_R <- solo_binary_R %>%
      group_by(marker) %>%
      summarise(
        x = sum(R, na.rm = TRUE),
        n = sum(pheno %in% c("S", "I", "R", "SDD")),
        p = x / n,
        se = sqrt(p * (1 - p) / n),
        ci.lower = max(0, p - 1.96 * se),
        ci.upper = min(1, p + 1.96 * se)
      ) %>%
      mutate(category = "R")
  } else {
    solo_stats_R <- tibble(
      marker = character(), x = double(), n = integer(),
      p = double(), se = double(), ci.lower = double(),
      ci.upper = double(), category = character()
    )
  }
  
  solo_stats_R <- solo_stats_R %>% arrange(p)
  if (reverse_order) {solo_stats_R <- solo_stats_R %>% arrange(desc(p))}

  if (icat & sum(!is.na(solo_binary$pheno)) > 0) {
    if ("I" %in% excludeRanges & "mic" %in% colnames(solo_binary)) {
      solo_binary_I <- solo_binary_norange
    } else {
      solo_binary_I <- solo_binary
    }
    solo_stats_I <- solo_binary_I %>%
      group_by(marker) %>%
      summarise(
        x = sum(I, na.rm = TRUE),
        n = sum(pheno %in% c("S", "I", "R", "SDD")),
        p = x / n,
        se = sqrt(p * (1 - p) / n),
        ci.lower = max(0, p - 1.96 * se),
        ci.upper = min(1, p + 1.96 * se)
      ) %>%
      mutate(category = "I")
  } else {
    solo_stats_I <- tibble(
      marker = character(), x = double(), n = integer(),
      p = double(), se = double(), ci.lower = double(),
      ci.upper = double(), category = character()
    )
  }

  if (sum(!is.na(solo_binary$ecoff)) > 0) {
    if ("NWT" %in% excludeRanges & "mic" %in% colnames(solo_binary)) {
      solo_binary_NWT <- solo_binary_norange
    } else {
      solo_binary_NWT <- solo_binary
    }
    solo_stats_NWT <- solo_binary_NWT %>%
      group_by(marker) %>%
      summarise(
        x = sum(NWT, na.rm = TRUE),
        n = sum(ecoff %in% c("WT", "NWT", "S", "I", "R", "SDD")),
        p = x / n,
        se = sqrt(p * (1 - p) / n),
        ci.lower = max(0, p - 1.96 * se),
        ci.upper = min(1, p + 1.96 * se)
      ) %>%
      mutate(category = "NWT")
  } else {
    solo_stats_NWT <- tibble(
      marker = character(), x = double(), n = integer(),
      p = double(), se = double(), ci.lower = double(),
      ci.upper = double(), category = character()
    )
  }

  solo_stats <- bind_rows(solo_stats_R, solo_stats_I) %>%
    bind_rows(solo_stats_NWT) %>%
    relocate(category, .before = x) %>%
    rename(ppv = p)
  
  

  # plots

  # markers with minimum number for either pheno or ecoff
  markers_to_plot <- solo_stats %>%
    filter(n >= min) %>%
    pull(marker) %>%
    unique()

  if (sum(!is.na(solo_binary$pheno)) > 0) { # if we have S/I/R call, plot it
    solo_pheno_plot <- solo_binary_R %>%
      filter(marker %in% markers_to_plot) %>%
      filter(!is.na(pheno)) %>%
      ggplot(aes(y = marker, fill = pheno)) +
      geom_bar(stat = "count", position = "fill") +
      scale_fill_sir(colours_SIR = colours_SIR) +
      geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(vjust = .5), size = 3) +
      scale_y_discrete(limits = markers_to_plot) +
      theme_light() +
      theme(
        axis.text.x = element_text(size = axis_label_size),
        axis.text.y = element_text(size = axis_label_size)
      ) +
      labs(y = "", x = "Proportion", fill = "Phenotype")
  } else if (sum(!is.na(solo_binary$ecoff)) > 0) { # if we only have ECOFF call, plot that instead
    solo_pheno_plot <- solo_binary_NWT %>%
      filter(marker %in% markers_to_plot) %>%
      filter(!is.na(ecoff)) %>%
      ggplot(aes(y = marker, fill = ecoff)) +
      geom_bar(stat = "count", position = "fill") +
      scale_fill_sir(colours_SIR = colours_SIR) +
      geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(vjust = .5), size = 3) +
      scale_y_discrete(limits = markers_to_plot) +
      theme_light() +
      theme(
        axis.text.x = element_text(size = axis_label_size),
        axis.text.y = element_text(size = axis_label_size)
      ) +
      labs(y = "", x = "Proportion", fill = "ECOFF")
  }

  labels <- solo_stats %>%
    filter(marker %in% markers_to_plot) %>%
    select(marker, category, n) %>%
    pivot_wider(id_cols = marker, names_from = category, values_from = n)

  if ("R" %in% colnames(labels) & "NWT" %in% colnames(labels)) {
    labels_vector <- paste0("(n=", labels$R, ",", labels$NWT, ")")
  } else if ("R" %in% colnames(labels)) {
    labels_vector <- paste0("(n=", labels$R, ")")
  } else if ("NWT" %in% colnames(labels)) {
    labels_vector <- paste0("(n=", labels$NWT, ")")
  } else {
    labels_vector <- rep("", nrow(labels))
  }

  names(labels_vector) <- labels$marker

  suppressWarnings(solo_stats_toplot <- solo_stats %>%
    filter(marker %in% markers_to_plot) %>%
    mutate(category = forcats::fct_relevel(category, "NWT", after = Inf)))

  ppv_plot <- solo_stats_toplot %>%
    ggplot(aes(y = marker, group = category, col = category)) +
    geom_vline(xintercept = 0.5, linetype = 2) +
    geom_linerange(aes(xmin = ci.lower, xmax = ci.upper), position = pd) +
    geom_point(aes(x = ppv), position = pd) +
    theme_bw() +
    scale_y_discrete(limits = markers_to_plot, labels = labels_vector[markers_to_plot], position = "right") +
    labs(y = "", x = "Solo PPV", col = "Category") +
    scale_colour_manual(values = colours_ppv) +
    theme(
      axis.text.x = element_text(size = axis_label_size),
      axis.text.y = element_text(size = axis_label_size)
    ) +
    xlim(0, 1)

  if (!is.null(drug_class_list)) {
    header <- paste("Solo markers for class:", paste0(drug_class_list, collapse = ", "))
  } else {
    header <- "Solo markers"
  }

  if (!is.null(antibiotic)) {
    subtitle <- paste("vs phenotype for drug:", antibiotic)
  } else {
    subtitle <- "vs drug phenotype"
  }

  combined_plot <- solo_pheno_plot + ggtitle(header, subtitle = subtitle) +
    ppv_plot +
    plot_layout(axes = "collect", guides = "collect")

  print(combined_plot)

  return(list(
    solo_stats = solo_stats, combined_plot = combined_plot,
    solo_binary = solo_binary, solo_binary_norange = solo_binary_norange,
    amr_binary = binary_matrix, plot_order = labels_vector[markers_to_plot]
  ))
}
