# Generate a Series of Plots for AMR Gene and Combination Analysis

This function generates a set of visualizations to analyze AMR gene
combinations, MIC values, and gene prevalence from an input geno-pheno
binary matrix. It creates several plots, including assay distributions,
phenotype breakdown and positive predictive values for each marker
combination. The
[amr_upset](https://AMRverse.github.io/AMRgen/reference/amr_upset.md)
and [ppv](https://AMRverse.github.io/AMRgen/reference/ppv.md) functions
can be used to generate standard data visualisations using the component
plots.

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

- assay:

  A character string indicating whether to plot MIC or disk diffusion
  data. Must be one of:

  - "mic": plot MIC data stored in column `mic`

  - "disk": plot disk diffusion data stored in column `disk`

- print_set_size:

  Logical indicating whether, if `plot_set_size` is TRUE, to print the
  number of strains with each marker combination on the plot. Default is
  FALSE.

- print_category_counts:

  Logical indicating whether, if `plot_category` is TRUE, to print the
  number of strains in each resistance category for each marker
  combination in the plot. Default is FALSE.

- SIR_col:

  A named vector of colours for the percentage bar plot. The names
  should be the phenotype categories (e.g., "R", "I", "S"), and the
  values should be valid color names or hexadecimal color codes. Default
  values are those used in the AMR package `scale_colour_sir()`.

- boxplot_col:

  Colour for lines of the box plots summarising the MIC distribution for
  each marker combination. Default is "grey".

- colours_ppv:

  A named vector of colours for the plot of PPV estimates. The names
  should be "R", "I" and "NWT", and the values should be valid color
  names or hexadecimal color codes.

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

- pd:

  Position dodge, i.e. spacing for the R/NWT values to be positioned
  above/below the line in the PPV plot. Default 'position_dodge(width =
  0.8)'.

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

combo_stats(binary_matrix, min_set_size = 3, order = "value", assay = "mic")
#> $summary
#> # A tibble: 103 × 19
#>    marker_list        marker_count     n combination_id   R.n   R.ppv R.ci_lower
#>    <chr>                     <dbl> <int> <chr>          <dbl>   <dbl>      <dbl>
#>  1 ""                            0  2590 0_0_0_0_0_0_0…    10 0.00386    0.00147
#>  2 "qnrB"                        1     1 0_0_0_0_0_0_0…     1 1          1      
#>  3 "parE_E460K, gyrA…            2     1 0_0_0_0_0_0_0…     1 1          1      
#>  4 "parE_D475E"                  1    61 0_0_0_0_0_0_0…     0 0          0      
#>  5 "qnrA1"                       1     2 0_0_0_0_0_0_0…     0 0          0      
#>  6 "gyrA_S83A"                   1     3 0_0_0_0_0_0_0…     0 0          0      
#>  7 "qnrB4"                       1     2 0_0_0_0_0_0_0…     2 1          1      
#>  8 "parE_I355T"                  1    24 0_0_0_0_0_0_0…     0 0          0      
#>  9 "marR_S3N"                    1    38 0_0_0_0_0_0_0…     4 0.105      0.00769
#> 10 "marR_S3N, parE_D…            2     4 0_0_0_0_0_0_0…     0 0          0      
#> # ℹ 93 more rows
#> # ℹ 12 more variables: R.ci_upper <dbl>, NWT.n <dbl>, NWT.ppv <dbl>,
#> #   NWT.ci_lower <dbl>, NWT.ci_upper <dbl>, median_excludeRangeValues <dbl>,
#> #   q25_excludeRangeValues <dbl>, q75_excludeRangeValues <dbl>,
#> #   n_excludeRangeValues <int>, median_ignoreRanges <dbl>,
#> #   q25_ignoreRanges <dbl>, q75_ignoreRanges <dbl>
#> 
#> $assay_plot

#> 
#> $setsize_plot

#> 
#> $marker_grid_plot

#> 
#> $marker_count_plot

#> 
#> $category_plot

#> 
#> $ppv_plot

#> 
```
