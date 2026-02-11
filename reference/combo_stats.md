# Generate a Series of Plots for AMR Gene and Combination Analysis

This function generates a set of visualizations to analyze AMR gene
combinations, MIC values, and gene prevalence from an input
genotype-phenotype binary matrix. It creates several plots, including
assay distributions, phenotype breakdown, and positive predictive values
for each marker combination. The
[`amr_upset()`](https://AMRverse.github.io/AMRgen/reference/amr_upset.md)
and [`ppv()`](https://AMRverse.github.io/AMRgen/reference/ppv.md)
functions can be used to generate standard data visualisations using the
component plots.

## Usage

``` r
combo_stats(
  binary_matrix,
  min_set_size = 2,
  order = "",
  assay = "mic",
  print_set_size = FALSE,
  print_category_counts = FALSE,
  SIR_col = c(S = "#3CAEA3", I = "#F6D55C", R = "#ED553B"),
  boxplot_col = "grey",
  colours_ppv = c(R = "maroon", NWT = "navy"),
  antibiotic = NULL,
  species = NULL,
  bp_site = NULL,
  guideline = "EUCAST 2025",
  bp_S = NULL,
  bp_R = NULL,
  ecoff_bp = NULL,
  pd = position_dodge(width = 0.8)
)
```

## Arguments

- binary_matrix:

  A data frame containing the original binary matrix output from the
  [`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md)
  function. Expected columns are an identifier (column 1, any name),
  `pheno` (class sir, with S/I/R categories to colour points), `mic`
  (class mic, with MIC values to plot), and other columns representing
  gene presence/absence (binary coded, i.e., 1 = present, 0 = absent).

- min_set_size:

  An integer specifying the minimum size for a gene set to be included
  in the analysis and plots. Default is 2. Only marker combinations with
  at least this number of occurrences are included in the plots.

- order:

  A character string indicating the order of the combinations on the
  x-axis. Options are:

  - `""` (default): decreasing frequency of combinations

  - `"genes"`: order by the number of genes in each combination

  - `"value"`: order by the median assay value (MIC or disk zone) for
    each combination

  - `"ppv"`: order by the PPV estimated for each combination

- assay:

  A character string indicating whether to plot MIC or disk diffusion
  data. Must be one of:

  - `"mic"` (default): plot MIC data stored in column `mic`

  - `"disk"`: plot disk diffusion data stored in column `disk`

  - `NULL`: don't plot or summarise assay data

- print_set_size:

  Logical indicating whether, if `plot_set_size=TRUE`, to print the
  number of strains with each marker combination on the plot. Default is
  `FALSE`.

- print_category_counts:

  Logical indicating whether, if `plot_category=TRUE`, to print the
  number of strains in each resistance category for each marker
  combination in the plot. Default is `FALSE`.

- SIR_col:

  A named vector of colours for the percentage bar plot. The names
  should be the phenotype categories (e.g., `"R"`, `"I"`, `"S"`), and
  the values should be valid color names or hexadecimal color codes.
  Default values are those used in the AMR package
  [`AMR::scale_colour_sir()`](https://amr-for-r.org/reference/plot.html).

- boxplot_col:

  Colour for lines of the box plots summarising the MIC distribution for
  each marker combination. Default is `"grey"`.

- colours_ppv:

  A named vector of colours for the plot of PPV estimates. The names
  should be `"R"`, `"I"` and `"NWT"`, and the values should be valid
  color names or hexadecimal color codes.

- antibiotic:

  Optional. Antibiotic name used to retrieve clinical breakpoints for
  annotation of the assay distribution plot.

- species:

  Optional. Species name used for breakpoint lookup.

- bp_site:

  Optional. Breakpoint site (e.g. "Non-meningitis") used when retrieving
  clinical breakpoints.

- guideline:

  Guideline used for breakpoint lookup. Default is `"EUCAST 2025"`.

- bp_S:

  (optional) S breakpoint to add to plot (numerical).

- bp_R:

  (optional) R breakpoint to add to plot (numerical).

- ecoff_bp:

  (optional) ECOFF breakpoint to add to plot (numerical).

- pd:

  A
  [`ggplot2::position_dodge()`](https://ggplot2.tidyverse.org/reference/position_dodge.html)
  object controlling horizontal spacing of points and confidence
  intervals in the PPV plot. Default is `position_dodge(width = 0.8)`.

## Value

A list containing the following elements:

- `summary`: A data frame summarizing each marker combination observed,
  including median MIC (and interquartile range), number of resistant
  isolates, and positive predictive value for resistance.

- `assay_plot`: MIC or disk diffusion distribution per combination.

- `setsize_plot`: Bar plot showing number of isolates per combination.

- `marker_grid_plot`: Dot plot showing gene combinations.

- `marker_count_plot`: Gene prevalence plot.

- `category_plot`: Stacked phenotype proportion plot.

- `ppv_plot`: Positive predictive value plot with confidence intervals.

## Details

Marker combinations are defined as unique patterns of gene
presence/absence across isolates. Only combinations observed in at least
`min_set_size` isolates are included in the analysis. Median and
interquartile range are calculated per combination. For MIC data,
medians may optionally exclude values expressed as ranges (e.g.
"\<=0.25").

## Examples

``` r
if (FALSE) { # \dontrun{
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

binary_matrix <- get_binary_matrix(
  geno_table = ecoli_geno,
  pheno_table = ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "pheno_clsi",
  keep_assay_values = TRUE,
  keep_assay_values_from = "mic"
)

combo_stats(binary_matrix, min_set_size = 3, order = "value", assay = "mic")
} # }
```
