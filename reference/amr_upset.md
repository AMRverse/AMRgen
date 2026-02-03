# Generate a Series of Plots for AMR Gene and Combination Analysis

This function generates a set of visualizations to analyze AMR gene
combinations, MIC values, and gene prevalence from a given binary
matrix. It creates several plots, including MIC distributions, a bar
plot for the percentage of strains per combination, a dot plot for gene
combinations in strains, and a plot for gene prevalence. It also outputs
a table summarizing the MIC distribution (median, lower, upper) and
number resistant, for each marker combination.

## Usage

``` r
amr_upset(
  binary_matrix,
  min_set_size = 2,
  order = "",
  plot_set_size = FALSE,
  plot_category = TRUE,
  print_category_counts = FALSE,
  print_set_size = FALSE,
  boxplot_col = "grey",
  assay = "mic",
  SIR_col = c(S = "#3CAEA3", SDD = "#8FD6C4", I = "#F6D55C", R = "#ED553B"),
  antibiotic = NULL,
  species = NULL,
  bp_site = NULL,
  guideline = "EUCAST 2025",
  bp_S = NULL,
  bp_R = NULL,
  ecoff_bp = NULL
)
```

## Arguments

- binary_matrix:

  A data frame containing the original binary matrix output from the
  `get_binary_matrix` function. Expected columns are an identifier
  (column 1, any name), `pheno` (class sir, with S/I/R categories to
  colour points), `mic` (class mic, with MIC values to plot), and other
  columns representing gene presence/absence (binary coded, i.e., 1 =
  present, 0 = absent).

- min_set_size:

  An integer specifying the minimum size for a gene set to be included
  in the analysis and plots. Default is 2. Only genes with at least this
  number of occurrences are included in the plots.

- order:

  A character string indicating the order of the combinations on the
  x-axis. Options are:

  - "" (default): decreasing frequency of combinations

  - "genes": order by the number of genes in each combination

  - "value": order by the median assay value (MIC or disk zone) for each
    combination.

- plot_set_size:

  Logical indicating whether to include a bar plot showing the set size
  (i.e., number of times each combination of markers is observed).
  Default is FALSE.

- plot_category:

  Logical indicating whether to include a stacked bar plot showing, for
  each marker combination, the proportion of samples with each phenotype
  classification (specified by the `pheno` column in the input file).
  Default is TRUE.

- print_category_counts:

  Logical indicating whether, if `plot_category` is TRUE, to print the
  number of strains in each resistance category for each marker
  combination in the plot. Default is FALSE.

- print_set_size:

  Logical indicating whether, if `plot_set_size` is TRUE, to print the
  number of strains with each marker combination on the plot. Default is
  FALSE.

- boxplot_col:

  Colour for lines of the box plots summarising the MIC distribution for
  each marker combination. Default is "grey".

- assay:

  A character string indicating whether to plot MIC or disk diffusion
  data. Must be one of:

  - "mic": plot MIC data stored in column `mic`

  - "disk": plot disk diffusion data stored in column `disk`

- SIR_col:

  A named vector of colours for the percentage bar plot. The names
  should be the phenotype categories (e.g., "R", "I", "S"), and the
  values should be valid color names or hexadecimal color codes. Default
  values are those used in the AMR package `scale_colour_sir()`.

- antibiotic:

  (optional) Name of antibiotic, so we can retrieve breakpoints to the
  assay value distribution plot.

- species:

  (optional) Name of species, so we can retrieve breakpoints to add to
  the assay value distribution plot.

- bp_site:

  (optional) Breakpoint site to retrieve (only relevant if also
  supplying `species` and `antibiotic` to retrieve breakpoints to plot,
  and not supplying breakpoints via `bp_S`, `bp_R`, `ecoff_bp`).

- guideline:

  (optional) Guideline to use when looking up breakpoints (default
  'EUCAST 2025')

- bp_S:

  (optional) S breakpoint to add to plot (numerical).

- bp_R:

  (optional) R breakpoint to add to plot (numerical).

- ecoff_bp:

  (optional) ECOFF breakpoint to add to plot (numerical).

## Value

A list containing the following elements:

- `plot`: A grid of plots displaying: (i) grid showing the marker
  combinations observed, MIC distribution per marker combination,
  frequency per marker and (optionally) phenotype classification and/or
  number of samples for each marker combination.

- `summary`: A data frame summarizing each marker combination observed,
  including median MIC (and interquartile range), number of resistant
  isolates, and positive predictive value for resistance.

## Examples

``` r
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
#>  Defining NWT in binary matrix using ecoff column provided: ecoff 

amr_upset(binary_matrix, min_set_size = 3, order = "value", assay = "mic")

#> $plot

#> 
#> $summary
#> # A tibble: 103 × 14
#>    marker_list             marker_count     n   R.ppv     R NWT.ppv   NWT
#>    <chr>                          <dbl> <int>   <dbl> <dbl>   <dbl> <dbl>
#>  1 ""                                 0  2590 0.00386    10  0.0197    51
#>  2 "qnrB"                             1     1 1           1  1          1
#>  3 "parE_E460K, gyrA_S83W"            2     1 1           1  1          1
#>  4 "parE_D475E"                       1    61 0           0  0          0
#>  5 "qnrA1"                            1     2 0           0  1          2
#>  6 "gyrA_S83A"                        1     3 0           0  0.667      2
#>  7 "qnrB4"                            1     2 1           2  1          2
#>  8 "parE_I355T"                       1    24 0           0  0          0
#>  9 "marR_S3N"                         1    38 0.105       4  0.263     10
#> 10 "marR_S3N, parE_D475E"             2     4 0           0  0.25       1
#> # ℹ 93 more rows
#> # ℹ 7 more variables: median_excludeRangeValues <dbl>,
#> #   q25_excludeRangeValues <dbl>, q75_excludeRangeValues <dbl>,
#> #   n_excludeRangeValues <int>, median_ignoreRanges <dbl>,
#> #   q25_ignoreRanges <dbl>, q75_ignoreRanges <dbl>
#> 
```
