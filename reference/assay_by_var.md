# Plot Assay Values Colored by a Variable

This function by deafault creates a stacked bar plot, where the x-axis
represents assay measurements (either MIC, Minimum Inhibitory
Concentration) or disk diffusion zones values), the y-axis indicates
their frequency, and the bars are colored by a variable indicated using
`colour_by` (by default, colours indicate whether the assay value is
expressed as a range or not). Plots can optionally be faceted on an
additional categorical variable. Optionally, the data can be plotted as
a grouped bar plot instead, with the assay measures shown on the y-axis,
grouped and coloured by the `colour_by` variable, by setting barplot =
`TRUE`. If breakpoints are provided, or species and drug are provided so
we can extract breakpoints, lines indicating the S/R breakpoints (solid
lines) and ECOFF (dashed line) will be added to the plot (and printed in
the subtitle).

## Usage

``` r
assay_by_var(
  pheno_table,
  pheno_drug = NULL,
  measure = "mic",
  colour_by = NULL,
  colours = NULL,
  facet_var = NULL,
  bp_site = NULL,
  bp_S = NULL,
  bp_R = NULL,
  bp_ecoff = NULL,
  species = NULL,
  guideline = "EUCAST 2025",
  bp_colours = c(S = "grey", R = "grey", E = "grey"),
  measure_axis_label = "Measurement",
  y_axis_label = "Count",
  colour_legend_label = NULL,
  plot_title = NULL,
  boxplot = FALSE,
  facet_nrow = NULL,
  facet_ncol = NULL
)
```

## Arguments

- pheno_table:

  Phenotype table in standard format as per import_pheno().

- pheno_drug:

  (optional) Name of a drug to filter the `drug` column, and to retrieve
  breakpoints for.

- measure:

  Name of the column with assay measurements to plot (default "mic").

- colour_by:

  (optional) Field name containing a variable to colour bars by (default
  NULL, which will colour each bar to indicate whether the value is
  expressed as a range or not).

- colours:

  (optional) Manual colour scale to use for bar plot. If NULL,
  `colour_by` variable is of class 'sir', bars will by default be
  coloured using standard SIR colours.

- facet_var:

  (optional) Column name containing a variable to facet on (default
  NULL).

- bp_site:

  (optional) Breakpoint site to retrieve (only relevant if also
  supplying `species` and `antibiotic` to retrieve breakpoints, and not
  supplying breakpoints via `bp_S`, `bp_R`, `ecoff`).

- bp_S:

  (optional) S breakpoint to plot.

- bp_R:

  (optional) R breakpoint to plot.

- bp_ecoff:

  (optional) ECOFF breakpoint to plot.

- species:

  (optional) Name of species, so we can retrieve breakpoints to print at
  the top of the plot to help interpret it.

- guideline:

  (optional) Guideline to use when looking up breakpoints (default
  'EUCAST 2025').

- bp_colours:

  (optional) Manual colour scale for breakpoint lines (default
  `c(S = "grey", R = "grey", E = "grey")`).

- measure_axis_label:

  (optional) String to label the measurement axis (x-axis for histogram,
  y-axis for boxplot, default `"Measurement"`).

- y_axis_label:

  (optional) String to label the y-axis in histogram plot (default
  `"Count"`).

- colour_legend_label:

  (optional) String to label the colour legend (default `NULL`, which
  results in plotting the variable name specified via the 'colour_by'
  parameter). Also used to label the x-axis if boxplot=`TRUE`.

- plot_title:

  (optional) String to title the plot (default indicates whether MIC or
  disk distribution is plotted, prefixed with the antibiotic name if
  provided, e.g. 'Ciprofloxacin MIC distribution')

- boxplot:

  (optional) If `TRUE`, plot the data as a grouped boxplot of assay
  measures, grouped and coloured by the `colour_by` variable. Summary
  statistics (median, geometric mean, and interquartile range of assay
  measures) are also computed, stratified by the `colour_by` and
  `facet_var` variables.

- facet_nrow:

  (optional) Number of rows for the facet grid (not used unless
  `facet_var` is provided).

- facet_ncol:

  (optional) Number of columns for the facet grid (not used unless
  `facet_var` is provided).

## Value

If boxplot=`FALSE`, the plot is returned as a single unnamed value. If
boxplot=`TRUE`, the plot is returned (\$plot) along with the summary
statistics (\$stats).

## Examples

``` r
# plot MIC distribution, highlighting values expressed as ranges
assay_by_var(
  pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
  measure = "mic"
)


# colour by SIR interpretation recorded in column 'pheno_clsi'
assay_by_var(
  pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
  measure = "mic", colour_by = "pheno_clsi"
)


# manually specify colours for the barplot
assay_by_var(
  pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
  measure = "mic", colour_by = "pheno_clsi",
  colours = c(S = "skyblue", I = "orange", R = "maroon")
)


# look up ECOFF and CLSI breakpoints and annotate these on the plot
assay_by_var(
  pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
  measure = "mic", colour_by = "pheno_clsi",
  species = "E. coli", guideline = "CLSI 2025"
)
#>   MIC breakpoints determined using AMR package: S <= 0.25 and R > 1


# facet by method
assay_by_var(
  pheno_table = ecoli_pheno, pheno_drug = "Ciprofloxacin",
  measure = "mic", colour_by = "pheno_clsi",
  species = "E. coli", guideline = "CLSI 2025",
  facet_var = "method"
)
#>   MIC breakpoints determined using AMR package: S <= 0.25 and R > 1

```
