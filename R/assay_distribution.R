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
#' Plot Assay Values Colored by a Variable
#'
#' This function by deafault creates a stacked bar plot, where the x-axis represents assay measurements (either MIC, Minimum Inhibitory Concentration) or disk diffusion zones values), the y-axis indicates their frequency, and the bars are colored by a variable indicated using `colour_by` (by default, colours indicate whether the assay value is expressed as a range or not). Plots can optionally be faceted on an additional categorical variable. Optionally, the data can be plotted as a grouped bar plot instead, with the assay measures shown on the y-axis, grouped and coloured by the `colour_by` variable, by setting barplot = `TRUE`. If breakpoints are provided, or species and drug are provided so we can extract breakpoints, lines indicating the S/R breakpoints (solid lines) and ECOFF (dashed line) will be added to the plot (and printed in the subtitle).
#' @param pheno_table Phenotype table in standard format as per import_pheno().
#' @param measure Name of the column with assay measurements to plot (default "mic").
#' @param colour_by (optional) Field name containing a variable to colour bars by (default NULL, which will colour each bar to indicate whether the value is expressed as a range or not).
#' @param colours (optional) Manual colour scale to use for bar plot. If NULL, `colour_by` variable is of class 'sir', bars will by default be coloured using standard SIR colours.
#' @param facet_var (optional) Column name containing a variable to facet on (default NULL).
#' @param pheno_drug (optional) Name of a drug to filter the `drug` column, and to retrieve breakpoints for.
#' @param species (optional) Name of species, so we can retrieve breakpoints to print at the top of the plot to help interpret it.
#' @param boxplot (optional) If `TRUE`, plot the data as a grouped boxplot of assay measures, grouped and coloured by the `colour_by` variable. Summary statistics (median, geometric mean, and interquartile range of assay measures) are also computed, stratified by the `colour_by` and `facet_var` variables.
#' @param bp_site (optional) Breakpoint site to retrieve (only relevant if also supplying `species` and `antibiotic` to retrieve breakpoints, and not supplying breakpoints via `bp_S`, `bp_R`, `ecoff`).
#' @param bp_S (optional) S breakpoint to plot.
#' @param bp_R (optional) R breakpoint to plot.
#' @param bp_ecoff (optional) ECOFF breakpoint to plot.
#' @param bp_colours (optional) Manual colour scale for breakpoint lines (default `c(S = "grey", R = "grey", E = "grey")`).
#' @param guideline (optional) Guideline to use when looking up breakpoints (default 'EUCAST 2025').
#' @param measure_axis_label (optional) String to label the measurement axis (x-axis for histogram, y-axis for boxplot, default `"Measurement"`).
#' @param y_axis_label (optional) String to label the y-axis in histogram plot (default `"Count"`).
#' @param colour_legend_label (optional) String to label the colour legend (default `NULL`, which results in plotting the variable name specified via the 'colour_by' parameter). Also used to label the x-axis if boxplot=`TRUE`.
#' @param plot_title (optional) String to title the plot (default indicates whether MIC or disk distribution is plotted, prefixed with the antibiotic name if provided, e.g. 'Ciprofloxacin MIC distribution')
#' @param facet_nrow (optional) Number of rows for the facet grid (not used unless `facet_var` is provided).
#' @param facet_ncol (optional) Number of columns for the facet grid (not used unless `facet_var` is provided).
#' @importFrom ggplot2 aes element_text facet_wrap geom_bar geom_vline ggplot labs theme scale_fill_manual sym geom_jitter
#' @importFrom methods is
#' @return If boxplot=`FALSE`, the plot is returned as a single unnamed value. If boxplot=`TRUE`, the plot is returned ($plot) along with the summary statistics ($stats).
#' @examples
#' # plot MIC distribution, highlighting values expressed as ranges
#' assay_by_var(
#'   pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
#'   measure = "mic"
#' )
#'
#' # colour by SIR interpretation recorded in column 'pheno_clsi'
#' assay_by_var(
#'   pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
#'   measure = "mic", colour_by = "pheno_clsi"
#' )
#'
#' # manually specify colours for the barplot
#' assay_by_var(
#'   pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
#'   measure = "mic", colour_by = "pheno_clsi",
#'   colours = c(S = "skyblue", I = "orange", R = "maroon")
#' )
#'
#' # look up ECOFF and CLSI breakpoints and annotate these on the plot
#' assay_by_var(
#'   pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
#'   measure = "mic", colour_by = "pheno_clsi",
#'   species = "E. coli", guideline = "CLSI 2025"
#' )
#'
#' # facet by method
#' assay_by_var(
#'   pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
#'   measure = "mic", colour_by = "pheno_clsi",
#'   species = "E. coli", guideline = "CLSI 2025",
#'   facet_var = "method"
#' )
#'
#' @export
assay_by_var <- function(pheno_table, pheno_drug = NULL, measure = "mic",
                         colour_by = NULL, colours = NULL, facet_var = NULL,
                         bp_site = NULL, bp_S = NULL, bp_R = NULL, bp_ecoff = NULL,
                         species = NULL, guideline = "EUCAST 2025",
                         bp_colours = c(S = "grey", R = "grey", E = "grey"),
                         measure_axis_label = "Measurement", y_axis_label = "Count",
                         colour_legend_label = NULL, plot_title = NULL,
                         boxplot = FALSE, facet_nrow = NULL, facet_ncol = NULL) {
  if (!is.null(pheno_drug)) {
    if ("drug" %in% colnames(pheno_table)) {
      pheno_table <- pheno_table %>% filter(drug == as.ab(pheno_drug))
      if (nrow(pheno_table) == 0) {
        stop("Drug '", pheno_drug, "' not found in drug column")
      }
    } else {
      warning("Column 'drug' not found in phenotype table, so can't filter to the specified pheno_drug.\nEnsure your input table is already filtered to the relevant drug.")
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
      warning("Colour variable '", colour_by, "' not found in input table")
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
    colours <- c(range = "maroon", value = "navy", `NA` = "grey")
  }

  # if species and pheno_drug are provided, but breakpoints aren't, check breakpoints to annotate plot
  if (measure %in% c("mic", "disk") & !is.null(species) & !is.null(pheno_drug)) {
    if (is.null(bp_S) | is.null(bp_R)) {
      breakpoints <- safe_execute(checkBreakpoints(species = species, guide = guideline, antibiotic = pheno_drug, bp_site = bp_site, assay = toupper(measure)))
      if (is.null(bp_R)) {
        bp_R <- breakpoints$breakpoint_R
      }
      if (is.null(bp_S)) {
        bp_S <- breakpoints$breakpoint_S
      }
    }
    if (is.null(bp_ecoff)) {
      bp_ecoff <- safe_execute(getBreakpoints(species = species, guide = "EUCAST 2025", antibiotic = pheno_drug, "ECOFF") %>%
        filter(method == toupper(measure)) %>% pull(breakpoint_S))
    }
  }

  pheno_data <- pheno_table %>% count(factor(!!sym(measure)), !!sym(colour_by))
  colnames(pheno_data)[1] <- measure

  # create subtitle reporting breakpoints
  if (!is.null(bp_S) | !is.null(bp_R) | !is.null(bp_ecoff)) {
    if (measure == "mic" & grepl("EUCAST", guideline)) {
      subtitle <- paste(guideline, "S <=", bp_S, "R>", bp_R, "ECOFF:", bp_ecoff)
    } else if (measure == "mic") {
      subtitle <- paste(guideline, "S <=", bp_S, "R>=", bp_R, "ECOFF:", bp_ecoff)
    } else if (measure == "disk" & grepl("EUCAST", guideline)) {
      subtitle <- paste(guideline, "S >=", bp_S, "R<", bp_R, "ECOFF:", bp_ecoff)
    } else if (measure == "disk") {
      subtitle <- paste(guideline, "S >=", bp_S, "R<=", bp_R, "ECOFF:", bp_ecoff)
    }
  } else {
    subtitle <- NULL
  }

  # get x coordinates for breakpoints after converting to factor
  if (!is.null(bp_S) | !is.null(bp_R) | !is.null(bp_ecoff)) {
    assay_order <- pheno_table %>% count(factor(!!sym(measure)))
    colnames(assay_order)[1] <- measure
    bp_S_hist <- c(1:nrow(assay_order))[assay_order[, 1] == bp_S]
    bp_R_hist <- c(1:nrow(assay_order))[assay_order[, 1] == bp_R]
    bp_ecoff_hist <- c(1:nrow(assay_order))[assay_order[, 1] == round(bp_ecoff, 2)]
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
    if (!is.null(pheno_drug)) {
      plot_title <- paste(pheno_drug, plot_title)
    }
  }

  stats <- NULL
  if (nrow(pheno_table) > 0) {
    if (!boxplot) {
      plot_all <- pheno_table %>%
        ggplot(aes(x = factor(!!sym(measure)))) +
        labs(
          x = measure_axis_label, y = y_axis_label,
          fill = colour_legend_label, subtitle = subtitle,
          title = plot_title
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

      if (is(pheno_table[[colour_by]], "sir")) {
        plot_all <- plot_all +
          geom_bar(aes(fill = !!sym(colour_by))) # don't treat as factor, will colour automatically by SIR
      } else { # treat as factor and apply manual colours if provided
        plot_all <- plot_all +
          geom_bar(aes(fill = factor(!!sym(colour_by))))
      }

      if (!is.null(colours)) {
        plot_all <- plot_all + scale_fill_manual(values = colours)
      }

      if (!is.null(facet_var)) {
        if (is.null(facet_ncol)) {
          facet_ncol <- 1
        }
        if (pheno_table %>% filter(!is.na(get(facet_var))) %>% nrow() > 0) {
          plot_all <- plot_all + facet_wrap(~ get(facet_var), ncol = facet_ncol, nrow = facet_nrow, scales = "free_y")
        }
      }
    } else { # grouped boxplot instead
      plot_all <- pheno_table %>%
        ggplot(aes(x = factor(!!sym(colour_by)), y = !!sym(measure))) +
        geom_boxplot() +
        geom_jitter(aes(col = factor(!!sym(colour_by)))) +
        labs(
          y = measure_axis_label, x = colour_legend_label,
          colour = colour_legend_label, subtitle = subtitle,
          title = plot_title
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

      if (!is.null(colours)) {
        plot_all <- plot_all + scale_colour_manual(values = colours)
      }

      if (!is.null(facet_var)) {
        if (pheno_table %>% filter(!is.na(get(facet_var))) %>% nrow() > 0) {
          plot_all <- plot_all + facet_wrap(~ get(facet_var), ncol = facet_ncol, nrow = facet_nrow)
        }
      }

      # calculate summary stats
      if (is.null(facet_var)) {
        stats <- pheno_table %>%
          group_by(!!sym(colour_by)) %>%
          summarise(
            n = n(),
            median = median(!!sym(measure), na.rm = TRUE),
            geom_mean = 2^(mean(log2(mic), na.rm = TRUE)),
            q25 = stats::quantile(!!sym(measure), 0.25, na.rm = TRUE),
            q75 = stats::quantile(!!sym(measure), 0.75, na.rm = TRUE)
          ) %>%
          ungroup()
        colnames(stats)[1] <- colour_by
      } else {
        stats <- pheno_table %>%
          group_by(!!sym(colour_by), !!sym(facet_var)) %>%
          summarise(
            n = n(),
            median = median(!!sym(measure), na.rm = TRUE),
            geom_mean = 2^(mean(log2(mic), na.rm = TRUE)),
            q25 = stats::quantile(!!sym(measure), 0.25, na.rm = TRUE),
            q75 = stats::quantile(!!sym(measure), 0.75, na.rm = TRUE)
          ) %>%
          ungroup()
        colnames(stats)[1] <- colour_by
        colnames(stats)[2] <- facet_var
      }
    }

    # add breakpoints to plot
    if (!is.null(bp_S)) {
      if (!boxplot) {
        plot_all <- plot_all + geom_vline(xintercept = bp_S_hist, colour = bp_colours["S"], linetype = 1)
      } else {
        plot_all <- plot_all + geom_hline(yintercept = bp_S, colour = bp_colours["S"], linetype = 1)
      }
    }
    if (!is.null(bp_R)) {
      if (!boxplot) {
        plot_all <- plot_all + geom_vline(xintercept = bp_R_hist, colour = bp_colours["R"], linetype = 1)
      } else {
        plot_all <- plot_all + geom_hline(yintercept = bp_R, colour = bp_colours["R"], linetype = 1)
      }
    }
    if (!is.null(bp_ecoff)) {
      if (!boxplot) {
        plot_all <- plot_all + geom_vline(xintercept = bp_ecoff_hist, colour = bp_colours["E"], linetype = 2)
      } else {
        plot_all <- plot_all + geom_hline(yintercept = bp_ecoff, colour = bp_colours["E"], linetype = 2)
      }
    }
  } else {
    plot_all <- NULL
  }

  plot_all

  if (!is.null(stats)) {
    return(list(plot = plot_all, stats = stats))
  } else {
    return(plot_all)
  }
}
