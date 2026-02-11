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
#
#' Generate a Stacked Bar Plot of Assay Values Colored by a Variable
#'
#' This function creates a stacked bar plot using `ggplot2`, where the x-axis represents MIC (Minimum Inhibitory Concentration) or disk values, the y-axis indicates their frequency, and the bars are colored by a variable (by default, colours indicate whether the assay value is expressed as a range or not). Plots can optionally be faceted on an additional categorical variable. If breakpoints are provided, or species and drug are provided so we can extract EUCAST breakpoints, vertical lines indicating the S/R breakpoints and ECOFF will be added to the plot.
#' @param pheno_table Phenotype table in standard format as per import_ast().
#' @param measure Name of the column with assay measurements to plot (default "mic").
#' @param colour_by (optional) Field name containing a variable to colour bars by (default NULL, which will colour each bar to indicate whether the value is expressed as a range or not).
#' @param bar_cols (optional) Manual colour scale to use for bar plot. If NULL, `colour_by` variable is of class 'sir', bars will by default be coloured using standard SIR colours.
#' @param facet_var (optional) Column name containing a variable to facet on (default NULL).
#' @param antibiotic (optional) Name of an antibiotic to filter the 'drug_agent' column, and to retrieve breakpoints for.
#' @param species (optional) Name of species, so we can retrieve breakpoints to print at the top of the plot to help interpret it.
#' @param bp_site (optional) Breakpoint site to retrieve (only relevant if also supplying `species` and `antibiotic` to retrieve breakpoints, and not supplying breakpoints via `bp_S`, `bp_R`, `ecoff`).
#' @param bp_S (optional) S breakpoint to plot.
#' @param bp_R (optional) R breakpoint to plot.
#' @param bp_ecoff (optional) ECOFF breakpoint to plot.
#' @param bp_cols (optional) Manual colour scale for breakpoint lines.
#' @param guideline (optional) Guideline to use when looking up breakpoints (default 'EUCAST 2025').
#' @param x_axis_label (optional) String to label the x-axis (default "Measurement").
#' @param y_axis_label (optional) String to label the y-axis (default "Count").
#' @param colour_legend_label (optional) String to label the barplot fill colour legend (default NULL, which results in plotting the variable name specified via the 'colour_by' parameter).
#' @param plot_title (optional) String to title the plot (default indicates whether MIC or disk distribution is plotted, prefixed with the antibiotic name if provided, e.g. 'Ciprofloxacin MIC distribution')
#' @importFrom ggplot2 aes element_text facet_wrap geom_bar geom_vline ggplot labs theme scale_fill_manual sym
#' @importFrom methods is
#' @return The stacked bar plot
#' @examples
#' # plot MIC distribution, highlighting values expressed as ranges
#' assay_by_var(
#'   pheno_table = ecoli_ast, antibiotic = "Ciprofloxacin",
#'   measure = "mic"
#' )
#'
#' # colour by SIR interpretation recorded in column 'pheno_clsi'
#' assay_by_var(
#'   pheno_table = ecoli_ast, antibiotic = "Ciprofloxacin",
#'   measure = "mic", colour_by = "pheno_clsi"
#' )
#'
#' # manually specify colours for the barplot
#' assay_by_var(
#'   pheno_table = ecoli_ast, antibiotic = "Ciprofloxacin",
#'   measure = "mic", colour_by = "pheno_clsi",
#'   bar_cols = c(S = "skyblue", I = "orange", R = "maroon")
#' )
#'
#' # look up ECOFF and CLSI breakpoints and annotate these on the plot
#' assay_by_var(
#'   pheno_table = ecoli_ast, antibiotic = "Ciprofloxacin",
#'   measure = "mic", colour_by = "pheno_clsi",
#'   species = "E. coli", guideline = "CLSI 2025"
#' )
#'
#' # facet by method
#' assay_by_var(
#'   pheno_table = ecoli_ast, antibiotic = "Ciprofloxacin",
#'   measure = "mic", colour_by = "pheno_clsi",
#'   species = "E. coli", guideline = "CLSI 2025",
#'   facet_var = "method"
#' )
#'
#' @export
assay_by_var <- function(pheno_table, antibiotic = NULL, measure = "mic",
                         colour_by = NULL, bar_cols = NULL, facet_var = NULL,
                         bp_site = NULL, bp_S = NULL, bp_R = NULL, bp_ecoff = NULL,
                         species = NULL, guideline = "EUCAST 2025",
                         bp_cols = c(S = "#3CAEA3", R = "#ED553B", E = "grey"),
                         x_axis_label = "Measurement", y_axis_label = "Count",
                         colour_legend_label = NULL, plot_title = NULL) {
  if (!is.null(antibiotic)) {
    if ("drug_agent" %in% pheno_table) {
      pheno_table <- pheno_table %>% filter(drug_agent == as.ab(antibiotic))
      if (nrow(pheno_table) == 0) {
        stop(paste0("Antibiotic '", antibiotic, "' not found in drug_agent column"))
      }
    }
  }

  if (measure %in% colnames(pheno_table)) {
    pheno_table <- pheno_table %>%
      filter(!is.na(get(measure))) %>%
      arrange(get(measure))
  } else {
    stop(paste0("No '", measure, "' column in input table"))
  }

  if (!is.null(facet_var)) {
    if (!(facet_var %in% colnames(pheno_table))) {
      stop(paste0("Facet variable '", facet_var, "' not found in input table"))
    }
  }

  if (!is.null(colour_by)) {
    if (!(colour_by %in% colnames(pheno_table))) {
      cat(paste0("WARNING: Colour variable '", colour_by, "' not found in input table\n"))
      colour_by <- NULL
    } else {
      pheno_table <- pheno_table %>%
        filter(!is.na(get(colour_by)))
    }
  }

  if (is.null(colour_by)) {
    pheno_table <- pheno_table %>%
      mutate(range = if_else(grepl("<", get(measure)), "range", "value"))
    colour_by <- "range"
    bar_cols <- c(range = "maroon", value = "navy", `NA` = "grey")
  }

  # if species and antibiotic are provided, but breakpoints aren't, check breakpoints to annotate plot
  if (measure %in% c("mic", "disk") & !is.null(species) & !is.null(antibiotic)) {
    if (is.null(bp_S) | is.null(bp_R)) {
      breakpoints <- safe_execute(checkBreakpoints(species = species, guide = guideline, antibiotic = antibiotic, bp_site = bp_site, assay = toupper(measure)))
      if (is.null(bp_R)) {
        bp_R <- breakpoints$breakpoint_R
      }
      if (is.null(bp_S)) {
        bp_S <- breakpoints$breakpoint_S
      }
    }
    if (is.null(bp_ecoff)) {
      bp_ecoff <- safe_execute(getBreakpoints(species = species, guide = "EUCAST 2025", antibiotic = antibiotic, "ECOFF") %>%
        filter(method == toupper(measure)) %>% pull(breakpoint_S))
    }
  }

  pheno_data <- pheno_table %>% count(factor(!!sym(measure)), !!sym(colour_by))
  colnames(pheno_data)[1] <- measure

  # create subtitle reporting breakpoints
  if (!is.null(bp_S) | !is.null(bp_R) | !is.null(bp_ecoff)) {
    if (measure == "mic" & grepl("EUCAST", guideline)) {
      subtitle <- paste("S <=", bp_S, "R>", bp_R, "ECOFF:", bp_ecoff)
    } else if (measure == "mic") {
      subtitle <- paste("S <=", bp_S, "R>=", bp_R, "ECOFF:", bp_ecoff)
    } else if (measure == "disk" & grepl("EUCAST", guideline)) {
      subtitle <- paste("S >=", bp_S, "R<", bp_R, "ECOFF:", bp_ecoff)
    } else if (measure == "disk") {
      subtitle <- paste("S >=", bp_S, "R<=", bp_R, "ECOFF:", bp_ecoff)
    }
  } else {
    subtitle <- NULL
  }

  # get x coordinates for breakpoints after converting to factor
  if (!is.null(bp_S) | !is.null(bp_R) | !is.null(bp_ecoff)) {
    assay_order <- pheno_table %>% count(factor(!!sym(measure)))
    colnames(assay_order)[1] <- measure
    bp_S <- c(1:nrow(assay_order))[assay_order[, 1] == bp_S]
    bp_R <- c(1:nrow(assay_order))[assay_order[, 1] == bp_R]
    bp_ecoff <- c(1:nrow(assay_order))[assay_order[, 1] == bp_ecoff]
  }

  # plot distribution per variable
  if (is.null(colour_legend_label)) {
    if (!is.null(colour_by)) {
      colour_legend_label <- colour_by
    }
  }
  if (is.null(plot_title)) {
    if (measure == "mic") {
      plot_title <- paste("MIC distribution")
    } else if (measure == "disk") {
      plot_title <- paste("Disk distribution")
    }
    if (!is.null(antibiotic)) {
      plot_title <- paste(antibiotic, plot_title)
    }
  }
  if (nrow(pheno_table) > 0) {
    plot_all <- pheno_table %>%
      ggplot(aes(x = factor(!!sym(measure)))) +
      labs(
        x = x_axis_label, y = y_axis_label,
        fill = colour_legend_label, subtitle = subtitle,
        title = plot_title
      ) +
      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
    if (is(pheno_table[[colour_by]], "sir")) {
      plot_all <- plot_all +
        geom_bar(aes(fill = !!sym(colour_by))) # don't treat as factor, will colour automatically by SIR
    } else { # treat as factor and apply manual colours if provided
      plot_all <- plot_all +
        geom_bar(aes(fill = factor(!!sym(colour_by))))
    }
    if (!is.null(bar_cols)) {
      plot_all <- plot_all + scale_fill_manual(values = bar_cols)
    }
    if (!is.null(facet_var)) {
      if (pheno_table %>% filter(!is.na(get(facet_var))) %>% nrow() > 0) {
        plot_all <- plot_all + facet_wrap(~ get(facet_var), ncol = 1, scales = "free_y")
      }
    }

    # add breakpoints to plot
    if (!is.null(bp_S)) {
      plot_all <- plot_all + geom_vline(xintercept = bp_S, color = bp_cols["S"], linetype = 2)
    }
    if (!is.null(bp_R)) {
      plot_all <- plot_all + geom_vline(xintercept = bp_R, color = bp_cols["R"], linetype = 2)
    }
    if (!is.null(bp_ecoff)) {
      plot_all <- plot_all + geom_vline(xintercept = bp_ecoff, color = bp_cols["E"], linetype = 2)
    }
  } else {
    plot_all <- NULL
  }

  return(plot_all)
}
