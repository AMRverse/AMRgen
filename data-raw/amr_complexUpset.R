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

#' Create an Upset plot for AMR data
#'
#' This function generates an Upset plot using a binary matrix of genetic determinants and antimicrobial resistance (AMR) data. Optionally, MIC or disk diffusion breakpoints can be overlaid on the plot.
#' @param binary_matrix A binary matrix where rows represent samples and columns represent genetic determinants or phenotypic resistance data.
#' @param min_set_size The minimum size of a set to be included in the plot. Default is 10.
#' @param mic_disk A character string specifying whether to use 'mic' or 'disk' for the y-axis. Must be either `"mic"` or `"disk"`. Default is `"mic"`.
#' @param remove_NAs Logical. If `TRUE`, removes rows with missing values in the selected `mic_disk` column. Default is `TRUE`.
#' @param gene_determinants A character vector specifying which genetic determinants to include. If `NULL`, all genes are used. Default is `NULL`.
#' @param colour_by A character string specifying the column used for colour mapping in the plot. Default is `"pheno"`.
#' @param plot_breakpoints Logical. If `TRUE`, overlays MIC or disk diffusion breakpoints on the plot. Default is `FALSE`.
#' @param organism A character string specifying the organism, used when plotting breakpoints.
#' @param break_guide A character string specifying the breakpoint guideline (e.g., `"EUCAST 2025"`). Default is `"EUCAST 2025"`.
#' @param break_type A character string specifying the breakpoint type (e.g., `"ECOFF"`). Default is `"ECOFF"`.
#' @param drug A character string specifying the antimicrobial drug to be analysed.
#' @param colour_values A named vector specifying colours for different resistance categories (`S`, `I`, `R`). Default is `c(S="#66c2a5", I="#fdae61", R="#d53e4f")`.
#' @param sort_intersections_by A character string passed to `ComplexUpset::upset()` to control sorting of intersections (e.g., `"degree"` or `"cardinality"`).
#' @param sort_intersections A character string specifying the order of intersections (`"ascending"` or `"descending"`).
#' @param show_intersect_size Logical. If `TRUE`, shows bar plot of intersection sizes. Default is `TRUE`.
#' @param intersect_counts Logical. If `TRUE`, counts are displayed in the intersection size bars. Default is `TRUE`.
#' @param heights Numeric vector specifying relative heights of the upset matrix vs annotation plots. If `NULL`, defaults based on `show_intersect_size`.
#' @importFrom AMR as.ab as.mic as.mo scale_y_mic
#' @importFrom ComplexUpset intersection_size upset
#' @importFrom dplyr across all_of filter mutate pull rename
#' @importFrom ggplot2 aes element_text geom_count geom_hline ggplot labs scale_colour_manual theme sym
#' @importFrom patchwork plot_layout
#' @return A `ggplot` object displaying the Upset plot.
#' @export
#' @examples
#' ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
#'
#' binary_matrix <- get_binary_matrix(
#'   geno_table = ecoli_geno,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno",
#'   keep_assay_values = TRUE,
#'   keep_assay_values_from = "mic"
#' )
#'
#' amr_complexUpset(binary_matrix)
amr_complexUpset <- function(binary_matrix,
                             min_set_size = 10,
                             mic_disk = "mic",
                             remove_NAs = TRUE, gene_determinants = NULL,
                             colour_by = "pheno",
                             plot_breakpoints = FALSE,
                             organism = NULL,
                             break_guide = "EUCAST 2025",
                             break_type = "ECOFF",
                             drug = NULL,
                             sort_intersections_by = "degree",
                             sort_intersections = "ascending",
                             show_intersect_size = TRUE,
                             intersect_counts = TRUE,
                             heights = NULL,
                             colour_values = c(S = "#66c2a5", I = "#fdae61", R = "#d53e4f")) {
  # mic_disk must be either 'mic' or 'disk'
  if (!mic_disk %in% c("mic", "disk")) {
    stop("mic_disk must be either 'mic' or 'disk'. Please select one of these values to continue.")
  }

  # Convert organism and antibiotic to mo and ab if they are provided
  if (!is.null(organism)) {
    organism <- as.mo(organism)
  }
  if (!is.null(drug)) {
    drug <- as.ab(drug)
  }

  # select gene determinants to plot
  if (is.null(gene_determinants)) {
    col <- colnames(binary_matrix)
    cols_to_remove <- c("mic", "disk", "R", "NWT", "pheno")
    genes <- col[-1]
    genes <- setdiff(genes, cols_to_remove)
  }
  # otherwise set to the vector provided
  else {
    genes <- gene_determinants
  }

  # convert gene columns to logical
  upset_data <- binary_matrix %>% mutate(across(all_of(genes), as.logical))

  # rename first column to 'sample' for the first column, to avoid issues with duplicate column names, if col1 is "id"
  if (colnames(upset_data)[1] == "id") {
    upset_data <- upset_data %>% rename(sample = id)
  }

  # if plot breakpoints is set, then get the relevant info (mo, ab, type and guideline MUST be set)
  if (plot_breakpoints) {
    # Check if all required arguments have values
    if (is.null(organism) || is.null(break_guide) || is.null(break_type) || is.null(drug)) {
      stop("When plot_breakpoints is TRUE, all of the following arguments must have values: organism, guideline, type, drug.")
    }
    # grab the method to filter by depending on whether we're using mic or disk
    if (mic_disk == "mic") {
      method_value <- "MIC"
    } else {
      method_value <- "DISK"
    }
    # extract the breakpoints
    breakpoint_table <- AMR::clinical_breakpoints %>%
      filter(ab == drug, mo == organism, guideline == break_guide, type == break_type, method == method_value)
    # if the table is empty, then return an error
    if (nrow(breakpoint_table) == 0) {
      warning("No breakpoints available for bug/drug combination. Will not plot breakpoint lines on upset, continuing...")
    }
    # otherwise grab the values, and create the hlines
    else {
      break_r <- breakpoint_table %>% pull(breakpoint_R)
      break_s <- breakpoint_table %>% pull(breakpoint_S)
      # convert to mic values if needed
      if (mic_disk == "mic") {
        break_r <- as.mic(break_r)
        break_s <- as.mic(break_s)
      }
      if (break_r == break_s) {
        warning("Breakpoint value is the same for R and S, only plotting one line")
        break_r_line <- geom_hline(yintercept = break_r, linetype = "dashed", alpha = 0.6, colour = "black")
        break_s_line <- NULL
      } else {
        break_r_line <- geom_hline(yintercept = as.mic(break_r), linetype = "dashed", alpha = 0.6, colour = "black")
        break_s_line <- geom_hline(yintercept = as.mic(break_s), linetype = "dashed", alpha = 0.6, colour = "black")
      }
    }
  } else {
    break_r_line <- NULL
    break_s_line <- NULL
  }

  # remove NAs if required
  if (remove_NAs) {
    upset_data <- upset_data %>% filter(!is.na(!!sym(mic_disk)))
  }

  # get the name of the y axis
  if (mic_disk == "mic") {
    y_axis_name <- "MIC (mg/L)"
    scale_y <- scale_y_mic()
  } else {
    y_axis_name <- "Disk (mm)"
    scale_y <- NULL
  }

  # check the colour palette matches the pheno values, otherwise create a new one
  col_vals <- names(table(upset_data[[colour_by]])[table(upset_data[[colour_by]]) > 0])
  if (!all(col_vals %in% names(colour_values))) {
    print(paste(
      "Warning: not all values of", colour_by, ":", paste(col_vals, collapse = ","),
      "are included in the colour vector `colour_values`:", paste(names(colour_values), collapse = ",")
    ))
    print("Defaulting to standard colours")
    # colour_values <-= rainbow(length(col_vals))
    colour_values <- stats::na.omit(colour_values[col_vals])
    to_add <- col_vals[!(col_vals %in% names(colour_values))]
    colour_values <- c(colour_values, grDevices::topo.colors(length(to_add)))
    names(colour_values) <- col_vals
  }

  if (show_intersect_size) {
    base_annotations <- list("Intersection size" = intersection_size(counts = intersect_counts))
    if (is.null(heights)) {
      heights <- c(2, 1)
    }
  } else {
    base_annotations <- list()
    if (is.null(heights)) {
      heights <- c(1, 1)
    }
  }

  # now plot
  plot <- upset(upset_data, genes,
    name = "genetic determinant", min_size = min_set_size, width_ratio = 0.1,
    sort_intersections_by = sort_intersections_by, sort_intersections = sort_intersections,
    base_annotations = base_annotations,
    annotations = list(
      y_axis_name = (
        ggplot(mapping = aes(y = !!sym(mic_disk), colour = as.factor(!!sym(colour_by)))) +
          geom_count() +
          break_r_line + # these values will be NULL if plot breakpoints isn't set
          break_s_line +
          scale_colour_manual(values = colour_values) +
          # scale_y + # set to scale_y_mic() if we have mic data
          labs(y = y_axis_name, colour = colour_by, size = "Number of\nisolates") +
          theme(legend.title = element_text(face = "bold"))
      )
    )
  ) + plot_layout(heights = heights) # relative heights of plotting areas

  return(plot)
}
