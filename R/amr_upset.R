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

#' Generate a Series of Plots for AMR Gene and Combination Analysis
#'
#' This function generates a set of visualizations to analyze AMR gene combinations, MIC values, and gene prevalence from an input genotype-phenotype binary matrix. It creates several plots, including assay distributions, phenotype breakdown, and positive predictive values for each marker combination. The [amr_upset()] and [ppv()] functions can be used to generate standard data visualisations using the component plots.
#' @param binary_matrix A data frame containing the original binary matrix output from the [get_binary_matrix()] function. Expected columns are an identifier (column 1, any name), `pheno` (class sir, with S/I/R categories to colour points), `mic` (class mic, with MIC values to plot), and other columns representing gene presence/absence (binary coded, i.e., 1 = present, 0 = absent).
#' @param min_set_size An integer specifying the minimum size for a gene set to be included in the analysis and plots. Default is 2. Only marker combinations with at least this number of occurrences are included in the plots.
#' @param order A character string indicating the order of the combinations on the x-axis. Options are:
#' - `""` (default): decreasing frequency of combinations
#' - `"genes"`: order by the number of genes in each combination
#' - `"value"`: order by the median assay value (MIC or disk zone) for each combination
#' - `"ppv"`: order by the PPV estimated for each combination
#' @param boxplot_col Colour for lines of the box plots summarising the MIC distribution for each marker combination. Default is `"grey"`.
#' @param SIR_col A named vector of colours for the percentage bar plot. The names should be the phenotype categories (e.g., `"R"`, `"I"`, `"S"`), and the values should be valid color names or hexadecimal color codes. Default values are those used in the AMR package [scale_colour_sir()].
#' @param assay A character string indicating whether to plot MIC or disk diffusion data. Must be one of:
#' - `"mic"` (default): plot MIC data stored in column `mic`
#' - `"disk"`: plot disk diffusion data stored in column `disk`
#' - `NULL`: don't plot or summarise assay data
#' @param print_set_size Logical indicating whether, if `plot_set_size=TRUE`, to print the number of strains with each marker combination on the plot. Default is `FALSE`.
#' @param print_category_counts Logical indicating whether, if `plot_category=TRUE`, to print the number of strains in each resistance category for each marker combination in the plot. Default is `FALSE`.
#' @param colours_ppv A named vector of colours for the plot of PPV estimates. The names should be `"R"`, `"I"` and `"NWT"`, and the values should be valid color names or hexadecimal color codes.
#' @param pd A `ggplot2::position_dodge()` object controlling horizontal spacing of points and confidence intervals in the PPV plot. Default is `position_dodge(width = 0.8)`.
#' @param bp_S (optional) S breakpoint to add to plot (numerical).
#' @param bp_R (optional) R breakpoint to add to plot (numerical).
#' @param ecoff_bp (optional) ECOFF breakpoint to add to plot (numerical).
#' @param antibiotic Optional. Antibiotic name used to retrieve clinical breakpoints for annotation of the assay distribution plot.
#' @param species Optional. Species name used for breakpoint lookup.
#' @param bp_site Optional. Breakpoint site (e.g. "Non-meningitis") used when retrieving clinical breakpoints.
#' @param guideline Guideline used for breakpoint lookup. Default is `"EUCAST 2025"`.
#' @importFrom AMR as.mic as.ab as.mo as.disk scale_color_sir scale_fill_sir scale_y_mic
#' @importFrom dplyr any_of arrange desc distinct filter group_by if_else left_join right_join mutate n pull relocate row_number select summarise ungroup
#' @importFrom forcats fct_rev fct_relevel
#' @importFrom ggplot2 aes after_stat coord_flip element_blank geom_bar geom_boxplot geom_col geom_point geom_segment geom_text ggplot labs position_fill scale_size_continuous scale_x_discrete scale_y_discrete scale_y_reverse theme theme_bw theme_light ylab ylim position_dodge geom_hline geom_linerange element_text scale_colour_manual
#' @importFrom patchwork plot_layout plot_spacer
#' @importFrom stats median quantile
#' @importFrom tidyr pivot_longer pivot_wider unite
#' @return A list containing the following elements:
#' - `summary`: A data frame summarizing each marker combination observed, including median MIC (and interquartile range), number of resistant isolates, and positive predictive value for resistance.
#' - `assay_plot`: MIC or disk diffusion distribution per combination.
#' - `setsize_plot`: Bar plot showing number of isolates per combination.
#' - `marker_grid_plot`: Dot plot showing gene combinations.
#' - `marker_count_plot`: Gene prevalence plot.
#' - `category_plot`: Stacked phenotype proportion plot.
#' - `ppv_plot`: Positive predictive value plot with confidence intervals.
#' @details
#' Marker combinations are defined as unique patterns of gene presence/absence
#' across isolates. Only combinations observed in at least `min_set_size`
#' isolates are included in the analysis. Median and interquartile range
#' are calculated per combination. For MIC data, medians may optionally
#' exclude values expressed as ranges (e.g. "<=0.25").
#' @export
#' @examples
#' \dontrun{
#' ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
#'
#' binary_matrix <- get_binary_matrix(
#'   geno_table = ecoli_geno,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno_clsi",
#'   keep_assay_values = TRUE,
#'   keep_assay_values_from = "mic"
#' )
#'
#' combo_matrix <- combo_stats(binary_matrix, min_set_size = 3, order = "value", assay = "mic")
#' }
combo_stats <- function(binary_matrix, min_set_size = 2, order = "",
                        assay = "mic",
                        print_set_size = FALSE,
                        print_category_counts = FALSE,
                        SIR_col = c(S = "#3CAEA3", I = "#F6D55C", R = "#ED553B"),
                        boxplot_col = "grey",
                        colours_ppv = c("R" = "maroon", "NWT" = "navy"),
                        antibiotic = NULL, species = NULL, bp_site = NULL,
                        guideline = "EUCAST 2025",
                        bp_S = NULL, bp_R = NULL, ecoff_bp = NULL,
                        pd = position_dodge(width = 0.8)) {
  if (sum(!is.na(binary_matrix$pheno)) == 0) {
    if (sum(!is.na(binary_matrix$ecoff)) > 0) {
      binary_matrix <- binary_matrix %>% mutate(pheno = ecoff)
      cat(" Warning: no values in pheno column, colouring upset plot by ecoff column\n")
    } else {
      stop(" Failed to make upset plot as no values in field pheno or ecoff")
    }
    colour_label <- "ECOFF\ncategory"
  } else {
    colour_label <- "Resistance\ncategory"
  }

  # remove rows with no call
  na_to_remove <- sum(is.na(binary_matrix$pheno))
  if (na_to_remove > 0) {
    cat(" Removing", na_to_remove, "rows with no phenotype call\n")
    binary_matrix <- binary_matrix %>% filter(!is.na(pheno))
  }

  # gene names
  genes <- binary_matrix %>%
    select(-any_of(c("id", "pheno", "ecoff", "R", "I", "NWT", "mic", "disk"))) %>%
    colnames()

  # check we have the expected assay column, with data
  if (!is.null(assay)) {
    if (!(assay %in% c("mic", "disk"))) {
      cat(paste("In valid `assay` value:", assay, "\n"))
      assay <- NULL
    } else if (!(assay %in% colnames(binary_matrix))) {
      stop(paste("input", deparse(substitute(binary_matrix)), "must have a column labelled ", assay))
    }
  }
  if (!is.null(assay)) {
    data_rows <- binary_matrix %>%
      filter(!is.na(get(assay))) %>%
      nrow()

    if (data_rows == 0) {
      stop(paste("input", deparse(substitute(binary_matrix)), "has no non-NA values in column ", assay))
    }

    # Add in a combination column and filter to samples with the required assay data (MIC or disk) only
    binary_matrix_wide <- get_combo_matrix(binary_matrix, assay = assay)
  } else {
    # Add in a combination column
    binary_matrix_wide <- get_combo_matrix(binary_matrix, assay = NULL)
  }

  # Make matrix longer
  binary_matrix <- binary_matrix_wide %>%
    pivot_longer(cols = genes[1]:genes[length(genes)], names_to = "genes") %>%
    mutate(genes = gsub("\\.\\.", ":", genes)) %>%
    mutate(genes = gsub("`", "", genes))

  ### Counts per combination, for bar plot - X axis = combination. Y axis = number of strains ###
  # This first to filter on combinations with enough data
  combination_freq <- binary_matrix_wide %>%
    group_by(combination_id) %>%
    summarise(n = n()) %>%
    mutate(perc = 100 * n / sum(n)) %>% # count number with each combination
    filter(n >= min_set_size) # count filter

  if (nrow(combination_freq) == 0) {
    stop(paste("No marker combinations pass the minimum frequency:", min_set_size))
  }
  # which have enough strains / data
  comb_enough_strains <- combination_freq %>% pull(combination_id)

  ### data for assay plot - dot plot. X axis = combination. Y axis = MIC or disk zone #####
  if (!is.null(assay)) {
    assay_plot <- binary_matrix %>%
      filter(combination_id %in% comb_enough_strains) %>%
      select(id:combination_id) %>%
      distinct() %>%
      group_by(combination_id, get(assay), pheno) %>%
      summarise(n = n()) # count how many at each assay value, keep pheno for colour
  } else {
    assay_plot <- NULL
  }

  ### Gene prevalence plot (amongst strains with marker combinations exceeding the minimum set size)
  gene.prev <- binary_matrix %>%
    filter(combination_id %in% comb_enough_strains) %>%
    group_by(genes) %>%
    summarise(gene.prev = sum(value))

  ### Set order of genes for y axis in combination dot plot
  ## (amongst strains with marker combinations exceeding the minimum set size, that will be included in plot)
  gene.order.desc <- gene.prev %>%
    arrange(desc(gene.prev)) %>%
    filter(gene.prev > 0) %>%
    pull(genes)

  # For gene prev plot
  gene.prev <- gene.prev %>%
    filter(genes %in% gene.order.desc) %>%
    mutate(genes = factor(genes, levels = gene.order.desc))

  # For combination dot plot
  binary_matrix <- binary_matrix %>%
    filter(combination_id %in% comb_enough_strains) %>%
    filter(genes %in% gene.order.desc) %>%
    mutate(genes = factor(genes, levels = gene.order.desc))

  ### Point plot - X axis = combination. Y axis = genes. Lines joining genes in same strain ###
  binary_matrix$point_size <- 2 * binary_matrix$value # want dot size to be larger than 1 => can make 2/3/4 etc

  # Only plot a line between points if more than one gene in a combination
  # get how many in each combination
  multi_genes_combination_id_all <- binary_matrix %>%
    select(combination_id, genes, value, point_size) %>%
    group_by(combination_id) %>%
    filter(value == 1) %>%
    mutate(u = length(unique(genes))) %>%
    filter(genes %in% gene.order.desc) %>%
    mutate(genes = factor(genes, levels = gene.order.desc, ordered = TRUE))

  # get only those with > 1
  multi_genes_combination_ids <- multi_genes_combination_id_all %>%
    filter(u > 1) %>%
    group_by(combination_id) %>%
    mutate(
      min = min(genes),
      max = max(genes)
    ) %>%
    ungroup()

  # summary table (ignore MIC values expressed as ranges, when calculating median/IQR)
  if (is.null(assay)) {
    summary <- binary_matrix_wide %>%
      group_by(combination_id) %>%
      summarise(n = n())
  } else if (assay == "mic") {
    summary <- binary_matrix_wide %>%
      group_by(combination_id) %>%
      summarise( # note these medians are not meaningful if all expressed as ranges
        median_ignoreRanges = median(mic, na.rm = TRUE),
        q25_ignoreRanges = stats::quantile(mic, 0.25, na.rm = TRUE),
        q75_ignoreRanges = stats::quantile(mic, 0.75, na.rm = TRUE),
        n = n()
      )
    summary <- binary_matrix_wide %>%
      filter(!grepl("<|>", as.character(mic))) %>%
      group_by(combination_id) %>%
      summarise(
        median_excludeRangeValues = median(mic, na.rm = TRUE),
        q25_excludeRangeValues = stats::quantile(mic, 0.25, na.rm = TRUE),
        q75_excludeRangeValues = stats::quantile(mic, 0.75, na.rm = TRUE),
        n_excludeRangeValues = n()
      ) %>%
      right_join(summary, by = "combination_id")
  } else if (assay == "disk") {
    summary <- binary_matrix_wide %>%
      group_by(combination_id) %>%
      summarise(
        median = median(as.double(as.disk(disk)), na.rm = TRUE),
        q25 = stats::quantile(as.double(as.disk(disk)), 0.25, na.rm = TRUE),
        q75 = stats::quantile(as.double(as.disk(disk)), 0.75, na.rm = TRUE),
        n = n(),
        n_exclRangeValues = sum(!grepl("<|>", as.character(mic)))
      )
  }
  if ("NWT" %in% colnames(binary_matrix_wide)) {
    summary <- binary_matrix_wide %>%
      group_by(combination_id) %>%
      summarise(NWT.n = sum(NWT, na.rm = TRUE), 
                NWT.denom = sum(!is.na(NWT))) %>%
      right_join(summary, by = "combination_id") %>%
      mutate(NWT.ppv = NWT.n / NWT.denom, .after = NWT.n) %>%
      mutate(NWT.se = sqrt(NWT.ppv * (1 - NWT.ppv) / NWT.denom)) %>%
      mutate(NWT.ci_upper = pmin(1, NWT.ppv + 1.96 * NWT.se), .after = NWT.ppv) %>%
      mutate(NWT.ci_lower = pmax(0, NWT.ppv - 1.96 * NWT.se), .after = NWT.ppv) %>%
      select(-NWT.se)
  }
  if ("I" %in% colnames(binary_matrix_wide)) {
    summary <- binary_matrix_wide %>%
      group_by(combination_id) %>%
      summarise(I.n = sum(I, na.rm = TRUE),
                I.denom = sum(!is.na(I))) %>%
      right_join(summary, by = "combination_id") %>%
      mutate(I.ppv = I.n / I.denom, .after = I.n) %>%
      mutate(I.se = sqrt(I.ppv * (1 - I.ppv) / I.denom)) %>%
      mutate(I.ci_upper = pmin(1, I.ppv + 1.96 * I.se), .after = I.ppv) %>%
      mutate(I.ci_lower = pmax(0, I.ppv - 1.96 * I.se), .after = I.ppv) %>%
      select(-I.se)
  }
  if ("R" %in% colnames(binary_matrix_wide)) {
    summary <- binary_matrix_wide %>%
      group_by(combination_id) %>%
      summarise(R.n = sum(R, na.rm = TRUE),
                R.denom = sum(!is.na(R))) %>%
      right_join(summary, by = "combination_id") %>%
      mutate(R.ppv = R.n / R.denom, .after = R.n) %>%
      mutate(R.se = sqrt(R.ppv * (1 - R.ppv) / R.denom)) %>%
      mutate(R.ci_upper = pmin(1, R.ppv + 1.96 * R.se), .after = R.ppv) %>%
      mutate(R.ci_lower = pmax(0, R.ppv - 1.96 * R.se), .after = R.ppv) %>%
      select(-R.se)
  }

  # get names for summary
  combination_names <- binary_matrix_wide %>%
    select(combination_id, any_of(genes)) %>%
    distinct()

  combination_names <- combination_names %>%
    mutate(marker_list = apply(., 1, function(row) {
      paste(names(combination_names)[-1][row[-1] == 1], collapse = ", ")
    }), .after = combination_id) %>%
    mutate(marker_count = rowSums(. == 1)) %>%
    select(combination_id, marker_list, marker_count)

  summary <- summary %>%
    left_join(combination_names, by = "combination_id") %>%
    mutate(marker_list = if_else(is.na(marker_list), "-", marker_list)) %>%
    mutate(marker_count = if_else(is.na(marker_count), 0, marker_count)) %>%
    relocate(marker_list, marker_count, n) %>%
    mutate(marker_list = gsub("\\.\\.", ":", marker_list)) %>%
    mutate(marker_list = gsub("`", "", marker_list))

  ### Set order of combination_id <- x axis
  # Default = decreasing frequency
  ordered_comb_order <- combination_freq %>%
    arrange(desc(perc)) %>%
    pull(combination_id)

  if (order == "ppv") {
    if ("R.ppv" %in% colnames(summary)) {
      ordered_comb_order <- summary %>%
        arrange(R.ppv) %>%
        pull(combination_id)
    } else if ("NWT.ppv" %in% colnames(summary)) {
      ordered_comb_order <- summary %>%
        arrange(R.ppv) %>%
        pull(combination_id)
    }
  }

  # Do by # genes in combination (only want each id once)
  else if (order == "genes") {
    ordered_comb_order <- multi_genes_combination_id_all %>%
      arrange(u) %>%
      filter(row_number() == 1) %>%
      pull(combination_id)
  }
  # Do by # median assay value in combination (only want each id once)
  else if (order == "value") {
    if (is.null(assay)) {
      cat("Order set to value, but `assay` not provided. Ordering by genes (set size) instead.\n")
    } else {
      if (assay == "mic") {
        # use mic class median (ie ignoring range indicators <>=) for the purpose of ordering columns
        mic_medians <- binary_matrix_wide %>%
          group_by(combination_id) %>%
          summarise(median = median(mic))
        ordered_comb_order <- mic_medians %>%
          arrange(median) %>%
          pull(combination_id)
      } else if (assay == "disk") {
        disk_medians <- binary_matrix_wide %>%
          group_by(combination_id) %>%
          summarise(median = median(as.double(disk)))
        ordered_comb_order <- disk_medians %>%
          arrange(-median) %>%
          pull(combination_id)
      }
    }
  }

  # apply the ordering
  combination_freq$combination_id <- factor(combination_freq$combination_id, levels = ordered_comb_order)
  binary_matrix$combination_id <- factor(binary_matrix$combination_id, levels = ordered_comb_order)
  summary$combination_id <- factor(summary$combination_id, levels = ordered_comb_order)
  if (!is.null(assay)) {
    assay_plot$combination_id <- factor(assay_plot$combination_id, levels = ordered_comb_order)
  }

  ##### Plots ###

  ### assay data plot (MIC/disk distribution)
  if (!is.null(assay)) {
    g1 <- ggplot(data = assay_plot, aes(x = combination_id, y = `get(assay)`)) +
      geom_boxplot(colour = boxplot_col) +
      geom_point(aes(size = n, colour = pheno), show.legend = TRUE) +
      theme_bw() +
      scale_size_continuous("Number of\nisolates") +
      scale_color_sir(
        colours_SIR = SIR_col,
        name = colour_label
      )

    if (assay == "mic") {
      g1 <- g1 +
        scale_y_mic() +
        ylab("MIC (mg/L)")
    } else {
      g1 <- g1 +
        ylab("Disk zone (mm)")
    }

    # if species and antibiotic are provided, but breakpoints aren't, check breakpoints to annotate plot
    if (!is.null(species) & !is.null(antibiotic) & (is.null(bp_S) | is.null(bp_R) | is.null(ecoff_bp))) {
      if (is.null(ecoff_bp)) {
        ecoff_bp <- safe_execute(getBreakpoints(species = as.mo(species), guide = "EUCAST 2025", antibiotic = as.ab(antibiotic), "ECOFF") %>% filter(method == toupper(assay)) %>% pull(breakpoint_S))
      }
      if (is.null(bp_S)) {
        bp_S <- safe_execute(unlist(checkBreakpoints(species = as.mo(species), guide = guideline, antibiotic = as.ab(antibiotic), bp_site = bp_site, assay = toupper(assay))[1]))
      }
      if (is.null(bp_R)) {
        bp_R <- safe_execute(unlist(checkBreakpoints(species = as.mo(species), guide = guideline, antibiotic = as.ab(antibiotic), bp_site = bp_site, assay = toupper(assay))[2]))
      }
    }

    if (!is.null(bp_S)) {
      g1 <- g1 + geom_hline(yintercept = bp_S)
    }
    if (!is.null(bp_R)) {
      g1 <- g1 + geom_hline(yintercept = bp_R)
    }
    if (!is.null(ecoff_bp)) {
      g1 <- g1 + geom_hline(yintercept = ecoff_bp, linetype = 2)
    }
  } else {
    g1 <- NULL
  }

  ### Bar plot - set size
  g2 <- ggplot(combination_freq, aes(x = combination_id, y = n)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_x_discrete("group")
  if (print_set_size) {
    g2 <- g2 + geom_text(aes(label = n), nudge_y = -.5, size = 3)
  }

  # pheno category stacked barplot
  category_plot <- binary_matrix %>%
    select(id, pheno, combination_id) %>%
    distinct() %>%
    ggplot(aes(x = combination_id, fill = pheno)) +
    geom_bar(stat = "count", position = "fill") +
    scale_fill_sir(colours_SIR = SIR_col) +
    theme_light() +
    labs(x = "", y = "Category", fill = "Resistance\ncategory")

  if (print_category_counts) {
    category_plot <- category_plot +
      geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(vjust = .5), size = 3)
  }

  ### Dot plot of combinations

  g3 <- binary_matrix <- binary_matrix %>%
    mutate(binary_comb = if_else(value > 0, 1, 0)) %>%
    ggplot(aes(x = combination_id, y = fct_rev(genes))) +
    geom_point(aes(size = binary_comb), show.legend = FALSE) +
    theme_bw() +
    scale_size_continuous(range = c(-1, 2)) +
    scale_y_discrete(name = "Marker") +
    geom_segment(
      data = multi_genes_combination_ids,
      aes(
        x = combination_id, xend = combination_id,
        y = min, yend = max, group = combination_id
      ),
      color = "black"
    )

  ### Plot gene prev / set size
  g4 <- gene.prev %>%
    ggplot(aes(x = fct_rev(genes), y = gene.prev)) +
    geom_col() +
    theme_bw() +
    coord_flip() +
    scale_y_reverse("")


  ppv_plot <- summary %>%
    rename(total = n) %>%
    filter(combination_id %in% comb_enough_strains) %>%
    pivot_longer(cols = starts_with(c("R", "I", "NWT")), names_to = c("category", "stat"), values_to = "value", names_pattern = "(.*)\\.(.*)") %>%
    pivot_wider(names_from = "stat", values_from = "value", id_cols = combination_id:category) %>%
    mutate(category = forcats::fct_relevel(category, "NWT", after = Inf)) %>%
    ggplot(aes(x = combination_id, group = category, col = category)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), position = pd) +
    geom_point(aes(y = ppv), position = pd) +
    theme_bw() +
    labs(x = "", y = "PPV", col = "Positive predictive\nvalue (PPV) for:") +
    scale_colour_manual(values = colours_ppv) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    ) +
    ylim(0, 1)


  return(list(
    summary = summary,
    assay_plot = g1,
    setsize_plot = g2,
    marker_grid_plot = g3,
    marker_count_plot = g4,
    category_plot = category_plot,
    ppv_plot = ppv_plot
  ))
}


#' Generate Upset Plot
#'
#' This function generates an upset plot showing summaries of phenotype results (assay distributions, phenotype category percentages, and/or predictive value for phenotype) for each combination of markers observed in the data.
#' @param binary_matrix A data frame containing the original binary matrix output from the [get_binary_matrix()] function. If not provided (or set to `NULL`), user must specify `geno_table`, `pheno_table`, `antibiotic`, `drug_class_list` and optionally `geno_sample_col`, `pheno_sample_col`, `sir_col`, `ecoff_col`, `marker_col` to pass to [get_binary_matrix()].
#' @param assay A character string indicating whether to plot MIC or disk diffusion data. Must be one of:
#' - `"mic"`: plot MIC data stored in column `mic`
#' - `"disk"`: plot disk diffusion data stored in column `disk`
#' @param min_set_size An integer specifying the minimum size for a gene set to be included in the analysis and plots. Default is 2. Only marker combinations with at least this number of occurrences are included in the plots.
#' @param order A character string indicating the order of the combinations on the x-axis. Options are:
#' - `""`: decreasing frequency of combinations
#' - `"genes"`: order by the number of genes in each combination
#' - `"value" (default)`: order by the median assay value (MIC or disk zone) for each combination
#' - `"ppv"`: order by the PPV estimated for each combination
#' @param geno_table (Required if `binary_matrix` not provided) A data frame containing genotype data, formatted with [import_amrfp()]. Only used if `binary_matrix` not provided.
#' @param pheno_table (Required if `binary_matrix` not provided) A data frame containing phenotype data, formatted with [import_ast()]. Only used if `binary_matrix` not provided.
#' @param antibiotic (Required if `binary_matrix` not provided) A character string specifying the antibiotic of interest to filter phenotype data (and optionally retrieve breakpoints). The value must match one of the entries in the `drug_agent` column of `pheno_table`. Only used if `binary_matrix` not provided or if breakpoints required.
#' @param drug_class_list (Required if `binary_matrix` not provided) A character vector of drug classes to filter genotype data for markers related to the specified antibiotic. Markers in `geno_table` will be filtered based on whether their `drug_class` matches any value in this list. Only used if `binary_matrix` not provided.
#' @param geno_sample_col A character string (optional) specifying the column name in `geno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers. Only used if `binary_matrix` not provided.
#' @param pheno_sample_col A character string (optional) specifying the column name in `pheno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers. Only used if `binary_matrix` not provided.
#' @param sir_col A character string specifying the column name in `pheno_table` that contains the resistance interpretation (SIR) data. The values should be `"S"`, `"I"`, `"R"` or otherwise interpretable by [AMR::as.sir()]. If not provided, the first column prefixed with "phenotype*" will be used if present, otherwise an error is thrown.  Only used if `binary_matrix` not provided.
#' @param ecoff_col A character string specifying the column name in `pheno_table` that contains resistance interpretations (SIR) made against the ECOFF rather than a clinical breakpoint. The values should be `"S"`, `"I"`, `"R"` or otherwise interpretable by [AMR::as.sir()]. Default `ecoff`. Set to `NULL` if not available.  Only used if `binary_matrix` not provided.
#' @param marker_col A character string specifying the column name in `geno_table` containing the marker identifiers. Default `"marker"`. Only used if `binary_matrix` not provided.
#' @param plot_set_size Logical indicating whether to include a bar plot showing the set size (i.e., number of times each combination of markers is observed). Default is `FALSE`.
#' @param print_set_size Logical indicating whether, if `plot_set_size=TRUE`, to print the number of strains with each marker combination on the plot. Default is `FALSE`.
#' @param plot_category Logical indicating whether to include a stacked bar plot showing, for each marker combination, the proportion of samples with each phenotype classification (specified by the `pheno` column in the input file). Default is `TRUE`.
#' @param print_category_counts Logical indicating whether, if `plot_category=TRUE`, to print the number of strains in each resistance category for each marker combination in the plot. Default is `FALSE`.
#' @param boxplot_col Colour for lines of the box plots summarising the MIC distribution for each marker combination. Default is `"grey"`.
#' @param SIR_col A named vector of colours for the percentage bar plot. The names should be the phenotype categories (e.g., `"R"`, `"I"`, `"S"`), and the values should be valid color names or hexadecimal color codes. Default values are those used in the AMR package [scale_colour_sir()].
#' @param bp_S (optional) S breakpoint to add to plot (numerical).
#' @param bp_R (optional) R breakpoint to add to plot (numerical).
#' @param ecoff_bp (optional) ECOFF breakpoint to add to plot (numerical).
#' @param antibiotic Optional. Antibiotic name used to retrieve clinical breakpoints for annotation of the assay distribution plot.
#' @param species Optional. Species name used for breakpoint lookup.
#' @param bp_site Optional. Breakpoint site (e.g. "Non-meningitis") used when retrieving clinical breakpoints.
#' @param guideline Guideline used for breakpoint lookup. Default is `"EUCAST 2025"`.
#' @importFrom ggplot2 theme theme_bw theme_light element_blank element_text labs scale_y_reverse
#' @importFrom patchwork plot_layout plot_spacer
#' @return A list containing the following elements:
#' - `plot`: A grid of plots displaying: (i) grid showing the marker combinations observed, MIC distribution per marker combination, frequency per marker and (optionally) phenotype classification and/or number of samples for each marker combination.
#' - `binary_matrix`: A copy of the genotype-phenotype binary matrix (either provided as input or generated by the function)
#' - `summary`: A data frame summarizing each marker combination observed, including median MIC (and interquartile range), number of resistant isolates, and positive predictive value for resistance.
#' @export
#' @examples
#' \dontrun{
#' ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
#'
#' # Generate binary matrix
#' binary_matrix <- get_binary_matrix(
#'   geno_table = ecoli_geno,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno_clsi",
#'   keep_assay_values = TRUE,
#'   keep_assay_values_from = "mic"
#' )
#'
#' # Run upset plot analysis using this binary_matrix
#' cip_mic_upset <- amr_upset(binary_matrix, assay = "mic")
#'
#' # Alternatively, generate binary matrix and run ppv() in one step
#' cip_mic_upset <- amr_upset(
#'   assay = "mic",
#'   geno_table = ecoli_geno,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno_clsi"
#' )
#' }
amr_upset <- function(binary_matrix = NULL, assay = "mic",
                      min_set_size = 2, order = "value",
                      geno_table, pheno_table, antibiotic = NULL, drug_class_list,
                      geno_sample_col = NULL, pheno_sample_col = NULL,
                      sir_col = NULL, ecoff_col = "ecoff",
                      marker_col = "marker",
                      plot_set_size = FALSE, plot_category = TRUE,
                      print_category_counts = FALSE, print_set_size = FALSE,
                      boxplot_col = "grey",
                      SIR_col = c(S = "#3CAEA3", I = "#F6D55C", R = "#ED553B"),
                      species = NULL, bp_site = NULL,
                      guideline = "EUCAST 2025",
                      bp_S = NULL, bp_R = NULL, ecoff_bp = NULL) {
  if (is.null(assay)) {
    stop("`assay` must be 'mic' or 'disk'.\nIf you don't want to plot assay data, function `ppv()` is more appropriate.")
  } else if (!(assay %in% c("mic", "disk"))) {
    stop("`assay` must be 'mic' or 'disk'.\nIf you don't want to plot assay data, function `ppv()` is more appropriate.")
  }

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

  combo_data <- combo_stats(
    binary_matrix = binary_matrix,
    min_set_size = min_set_size,
    order = order,
    assay = assay,
    boxplot_col = boxplot_col,
    print_set_size = print_set_size,
    print_category_counts = print_category_counts,
    SIR_col = SIR_col,
    antibiotic = antibiotic,
    species = species,
    bp_site = bp_site,
    guideline = guideline,
    bp_S = bp_S,
    bp_R = bp_R,
    ecoff_bp = ecoff_bp
  )

  # assemble plot
  final_plot <- plot_spacer() +
    combo_data$assay_plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    plot_layout(ncol = 2, widths = c(1, 4), guides = "collect")

  if (plot_category) {
    final_plot <- final_plot + plot_spacer() +
      combo_data$category_plot +
      theme(legend.position = "none") + # legend already done for assay plot
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
  }

  final_plot <- final_plot +
    combo_data$marker_count_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
    combo_data$marker_grid_plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
  if (plot_set_size) {
    final_plot <- final_plot + plot_spacer() +
      combo_data$setsize_plot + scale_y_reverse("Set size") +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
  }

  # set relative plotting heights
  if (plot_category & plot_set_size) {
    final_plot <- final_plot + plot_layout(heights = c(2, 1, 2, 1))
  }
  if (plot_category & !plot_set_size) {
    final_plot <- final_plot + plot_layout(heights = c(2, 1, 2))
  }
  if (!plot_category & plot_set_size) {
    final_plot <- final_plot + plot_layout(heights = c(2, 2, 1))
  }

  print(final_plot)

  return(list(plot = final_plot, binary_matrix = binary_matrix, summary = combo_data$summary))
}


#' Generate PPV Plot
#'
#' This function generates a plot showing summaries of phenotype results (assay distributions, phenotype category percentages, and/or predictive value for phenotype) for each combination of markers observed in the data.
#' @param binary_matrix A data frame containing the original binary matrix output from the [get_binary_matrix()] function. If not provided (or set to `NULL`), user must specify `geno_table`, `pheno_table`, `antibiotic`, `drug_class_list` and optionally `geno_sample_col`, `pheno_sample_col`, `sir_col`, `ecoff_col`, `marker_col` to pass to [get_binary_matrix()].
#' @param min_set_size An integer specifying the minimum size for a gene set to be included in the analysis and plots. Default is 2. Only marker combinations with at least this number of occurrences are included in the plots.
#' @param order A character string indicating the order of the combinations on the x-axis. Options are:
#' - `""`: decreasing frequency of combinations
#' - `"genes"`: order by the number of genes in each combination
#' - `"value"`: order by the median assay value (MIC or disk zone) for each combination.
#' - `"ppv"` (default): order by the PPV estimated for each combination
#' @param geno_table (Required if `binary_matrix` not provided) A data frame containing genotype data, formatted with [import_amrfp()]. Only used if `binary_matrix` not provided.
#' @param pheno_table (Required if `binary_matrix` not provided) A data frame containing phenotype data, formatted with [import_ast()]. Only used if `binary_matrix` not provided.
#' @param antibiotic (Required if `binary_matrix` not provided) A character string specifying the antibiotic of interest to filter phenotype data (and optionally retrieve breakpoints). The value must match one of the entries in the `drug_agent` column of `pheno_table`. Only used if `binary_matrix` not provided or if breakpoints required.
#' @param drug_class_list (Required if `binary_matrix` not provided) A character vector of drug classes to filter genotype data for markers related to the specified antibiotic. Markers in `geno_table` will be filtered based on whether their `drug_class` matches any value in this list. Only used if `binary_matrix` not provided.
#' @param geno_sample_col A character string (optional) specifying the column name in `geno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers. Only used if `binary_matrix` not provided.
#' @param pheno_sample_col A character string (optional) specifying the column name in `pheno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers. Only used if `binary_matrix` not provided.
#' @param sir_col A character string specifying the column name in `pheno_table` that contains the resistance interpretation (SIR) data. The values should be `"S"`, `"I"`, `"R"` or otherwise interpretable by [AMR::as.sir()]. If not provided, the first column prefixed with "phenotype*" will be used if present, otherwise an error is thrown.  Only used if `binary_matrix` not provided.
#' @param ecoff_col A character string specifying the column name in `pheno_table` that contains resistance interpretations (SIR) made against the ECOFF rather than a clinical breakpoint. The values should be `"S"`, `"I"`, `"R"` or otherwise interpretable by [AMR::as.sir()]. Default `ecoff`. Set to `NULL` if not available.  Only used if `binary_matrix` not provided.
#' @param marker_col A character string specifying the column name in `geno_table` containing the marker identifiers. Default `"marker"`. Only used if `binary_matrix` not provided.
#' @param colours_ppv A named vector of colours for the plot of PPV estimates. The names should be `"R"`, `"I"` and `"NWT"`, and the values should be valid color names or hexadecimal color codes.
#' @param SIR_col A named vector of colours for the percentage bar plot. The names should be the phenotype categories (e.g., `"R"`, `"I"`, `"S"`), and the values should be valid color names or hexadecimal color codes. Default values are those used in the AMR package [scale_colour_sir()].
#' @param upset_grid Logical indicating whether to show marker combinations as an upset plot-style grid (default `FALSE`, so that each row is instead labelled with a printed list of markers).
#' @param marker_label_space Relative width of plotting area to provide to the marker list/grid. (Default `NULL`, which results in a default value of 3 when `upset_grid=FALSE` and 1 otherwise).
#' @param plot_category Logical indicating whether to include a stacked bar plot showing, for each marker combination, the proportion of samples with each phenotype classification (specified by the `pheno` column in the input file). Default is `TRUE`.
#' @param print_category_counts Logical indicating whether, if `plot_category=TRUE`, to print the number of strains in each resistance category for each marker combination in the plot. Default is `FALSE`.
#' @param plot_ppv Logical indicating whether to plot the estimates for positive predictive value, for each marker combination (default `TRUE`).
#' @param plot_assay Logical indicating whether to plot the distribution of MIC/disk assay values, for each marker combination (default `FALSE`). Note you must also indicate which assay column to plot (`"mic"` or `"disk"`) via `assay`.
#' @param assay A character string indicating whether to plot MIC or disk diffusion data. Must be one of:
#' - `NULL`: (default) if no assay data is to be plotted
#' - `"mic"`: plot MIC data stored in column `mic`
#' - `"disk"`: plot disk diffusion data stored in column `disk`
#' @param boxplot_col Colour for lines of the box plots summarising the MIC distribution for each marker combination. Default is `"grey"`.
#' @param species Optional. Species name used to retrieve clinical breakpoints for annotation of the assay distribution plot.
#' @param bp_site Optional. Breakpoint site (e.g. "Non-meningitis") used when retrieving clinical breakpoints.
#' @param guideline Guideline used for breakpoint lookup. Default is `"EUCAST 2025"`.
#' @param bp_S (optional) S breakpoint to add to assay distribution plot (numerical).
#' @param bp_R (optional) R breakpoint to add to assay distribution plot (numerical).
#' @param ecoff_bp (optional) ECOFF breakpoint to add to assay distribution plot (numerical).
#' @param pd A `ggplot2::position_dodge()` object controlling horizontal spacing of points and confidence intervals in the PPV plot. Default is `position_dodge(width = 0.8)`.
#' @importFrom dplyr distinct filter group_by if_else left_join mutate n pull relocate row_number select summarise ungroup starts_with
#' @importFrom ggplot2 aes after_stat coord_flip element_blank geom_bar geom_boxplot geom_col geom_point geom_segment geom_text ggplot labs position_fill scale_size_continuous scale_x_discrete scale_y_discrete scale_y_reverse theme theme_bw theme_light ylab scale_x_continuous theme_void guides margin coord_cartesian
#' @importFrom patchwork plot_layout plot_spacer
#' @importFrom stats median
#' @importFrom tidyr pivot_longer unite
#' @return A list containing the following elements:
#' - `plot`: A grid of the requested plots
#' - `binary_matrix`: A copy of the genotype-phenotype binary matrix (either provided as input or generated by the function)
#' - `summary`: A data frame summarizing each marker combination observed, including number of resistant isolates, positive predictive values, and median assay values (and interquartile range) where relevant.
#' @export
#' @examples
#' \dontrun{
#' ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
#'
#' # Generate binary matrix
#' binary_matrix <- get_binary_matrix(
#'   geno_table = ecoli_geno,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno_clsi",
#'   keep_assay_values = TRUE,
#'   keep_assay_values_from = "mic"
#' )
#'
#' # Run ppv analysis using this binary_matrix
#' ppv <- ppv(binary_matrix)
#'
#' # Alternatively, generate binary matrix and run ppv() in one step
#' ppv <- ppv(
#'   geno_table = ecoli_geno,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno_clsi"
#' )
#' }
ppv <- function(binary_matrix = NULL,
                min_set_size = 2,
                order = "ppv",
                geno_table, pheno_table, antibiotic = NULL, drug_class_list,
                geno_sample_col = NULL, pheno_sample_col = NULL,
                sir_col = NULL, ecoff_col = "ecoff",
                marker_col = "marker",
                colours_ppv = c("R" = "maroon", "NWT" = "navy"),
                SIR_col = c(S = "#3CAEA3", I = "#F6D55C", R = "#ED553B"),
                upset_grid = FALSE,
                marker_label_space = NULL,
                plot_category = TRUE,
                print_category_counts = TRUE,
                plot_ppv = TRUE,
                plot_assay = FALSE,
                assay = NULL,
                boxplot_col = "grey", species = NULL, bp_site = NULL,
                guideline = "EUCAST 2025",
                bp_S = NULL, bp_R = NULL, ecoff_bp = NULL,
                pd = position_dodge(width = 0.8)) {
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

  combo_data <- combo_stats(
    binary_matrix = binary_matrix,
    min_set_size = min_set_size,
    order = order,
    assay = assay,
    colours_ppv = colours_ppv,
    boxplot_col = boxplot_col,
    print_category_counts = print_category_counts,
    SIR_col = SIR_col,
    antibiotic = antibiotic,
    species = species,
    bp_site = bp_site,
    guideline = guideline,
    bp_S = bp_S,
    bp_R = bp_R,
    ecoff_bp = ecoff_bp,
    pd = pd
  )

  # assemble plot
  if (upset_grid) {
    final_plot <- combo_data$marker_grid_plot + coord_flip() +
      scale_y_discrete(position = "right") +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0)
      ) +
      ylab("")
  } else {
    final_plot <- combo_data$summary %>%
      filter(n >= min_set_size) %>%
      mutate(marker_list = if_else(marker_list == "", "(none)", marker_list)) %>%
      ggplot(aes(y = combination_id, label = marker_list)) +
      geom_text(aes(x = Inf), hjust = 1, size = 3) +
      scale_x_continuous(limits = c(1, 1.5), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(
        plot.margin = margin(5.5, 1, 5.5, 0)
      )
  }

  if (plot_category) {
    final_plot <- final_plot +
      combo_data$category_plot +
      coord_flip() +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  if (plot_ppv) {
    final_plot <- final_plot +
      combo_data$ppv_plot +
      coord_flip() +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  if (plot_assay) {
    final_plot <- final_plot +
      combo_data$assay_plot +
      coord_flip() +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    if (plot_category) {
      final_plot <- final_plot + guides(color = "none") # already have SIR legend
    }
  }

  # add counts for each row to the right hand side
  
  if ("R.denom" %in% colnames(combo_data$summary) & "NWT.denom" %in% colnames(combo_data$summary)) {
    row_counts <- combo_data$summary %>%
      mutate(count_label = paste0("(n=", R.denom, ", ", NWT.denom, ")"))
  } else if ("R.denom" %in% colnames(combo_data$summary)) {
    row_counts <- combo_data$summary %>%
      mutate(count_label = paste0("(n=", R.denom, ")"))
  } else if ("NWT.denom" %in% colnames(combo_data$summary)) {
    row_counts <- combo_data$summary %>%
      mutate(count_label = paste0("(n=", NWT.denom, ")"))
  } 
      
  row_counts <- row_counts %>%
    filter(n >= min_set_size) %>%
    ggplot(aes(y = combination_id, x = 1, label = count_label)) +
    geom_text(hjust = 0, size = 3) +
    scale_x_continuous(limits = c(1, 1.5), expand = c(0, 0)) +
    theme_void() +
    theme(
      plot.margin = margin(5.5, 5.5, 5.5, 0)
    )

  if (is.null(marker_label_space)) {
    if (upset_grid) {
      marker_label_space <- 1
    } else {
      marker_label_space <- 3
    }
  }
  final_plot <- final_plot +
    row_counts +
    plot_layout(guides = "collect", nrow = 1, widths = c(marker_label_space, rep(1, sum(plot_category, plot_ppv, plot_assay)), 0.3))

  print(final_plot)

  return(list(plot = final_plot, binary_matrix = binary_matrix, summary = combo_data$summary))
}
