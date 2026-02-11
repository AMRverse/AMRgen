# Generate PPV Plot

This function generates a plot showing summaries of phenotype results
(assay distributions, phenotype category percentages, and/or predictive
value for phenotype) for each combination of markers observed in the
data.

## Usage

``` r
ppv(
  binary_matrix = NULL,
  min_set_size = 2,
  order = "ppv",
  geno_table,
  pheno_table,
  antibiotic = NULL,
  drug_class_list,
  geno_sample_col = NULL,
  pheno_sample_col = NULL,
  sir_col = NULL,
  ecoff_col = "ecoff",
  marker_col = "marker",
  colours_ppv = c(R = "maroon", NWT = "navy"),
  SIR_col = c(S = "#3CAEA3", I = "#F6D55C", R = "#ED553B"),
  upset_grid = FALSE,
  marker_label_space = NULL,
  plot_category = TRUE,
  print_category_counts = TRUE,
  plot_ppv = TRUE,
  plot_assay = FALSE,
  assay = NULL,
  boxplot_col = "grey",
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
  function. If not provided (or set to `NULL`), user must specify
  `geno_table`, `pheno_table`, `antibiotic`, `drug_class_list` and
  optionally `geno_sample_col`, `pheno_sample_col`, `sir_col`,
  `ecoff_col`, `marker_col` to pass to
  [`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md).

- min_set_size:

  An integer specifying the minimum size for a gene set to be included
  in the analysis and plots. Default is 2. Only marker combinations with
  at least this number of occurrences are included in the plots.

- order:

  A character string indicating the order of the combinations on the
  x-axis. Options are:

  - `""`: decreasing frequency of combinations

  - `"genes"`: order by the number of genes in each combination

  - `"value"`: order by the median assay value (MIC or disk zone) for
    each combination.

  - `"ppv"` (default): order by the PPV estimated for each combination

- geno_table:

  (Required if `binary_matrix` not provided) A data frame containing
  genotype data, formatted with
  [`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md).
  Only used if `binary_matrix` not provided.

- pheno_table:

  (Required if `binary_matrix` not provided) A data frame containing
  phenotype data, formatted with
  [`import_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ast.md).
  Only used if `binary_matrix` not provided.

- antibiotic:

  (Required if `binary_matrix` not provided) A character string
  specifying the antibiotic of interest to filter phenotype data (and
  optionally retrieve breakpoints). The value must match one of the
  entries in the `drug_agent` column of `pheno_table`. Only used if
  `binary_matrix` not provided or if breakpoints required.

- drug_class_list:

  (Required if `binary_matrix` not provided) A character vector of drug
  classes to filter genotype data for markers related to the specified
  antibiotic. Markers in `geno_table` will be filtered based on whether
  their `drug_class` matches any value in this list. Only used if
  `binary_matrix` not provided.

- geno_sample_col:

  A character string (optional) specifying the column name in
  `geno_table` containing sample identifiers. Defaults to `NULL`, in
  which case it is assumed the first column contains identifiers. Only
  used if `binary_matrix` not provided.

- pheno_sample_col:

  A character string (optional) specifying the column name in
  `pheno_table` containing sample identifiers. Defaults to `NULL`, in
  which case it is assumed the first column contains identifiers. Only
  used if `binary_matrix` not provided.

- sir_col:

  A character string specifying the column name in `pheno_table` that
  contains the resistance interpretation (SIR) data. The values should
  be `"S"`, `"I"`, `"R"` or otherwise interpretable by
  [`AMR::as.sir()`](https://amr-for-r.org/reference/as.sir.html). If not
  provided, the first column prefixed with "phenotype\*" will be used if
  present, otherwise an error is thrown. Only used if `binary_matrix`
  not provided.

- ecoff_col:

  A character string specifying the column name in `pheno_table` that
  contains resistance interpretations (SIR) made against the ECOFF
  rather than a clinical breakpoint. The values should be `"S"`, `"I"`,
  `"R"` or otherwise interpretable by
  [`AMR::as.sir()`](https://amr-for-r.org/reference/as.sir.html).
  Default `ecoff`. Set to `NULL` if not available. Only used if
  `binary_matrix` not provided.

- marker_col:

  A character string specifying the column name in `geno_table`
  containing the marker identifiers. Default `"marker"`. Only used if
  `binary_matrix` not provided.

- colours_ppv:

  A named vector of colours for the plot of PPV estimates. The names
  should be `"R"`, `"I"` and `"NWT"`, and the values should be valid
  color names or hexadecimal color codes.

- SIR_col:

  A named vector of colours for the percentage bar plot. The names
  should be the phenotype categories (e.g., `"R"`, `"I"`, `"S"`), and
  the values should be valid color names or hexadecimal color codes.
  Default values are those used in the AMR package
  [`AMR::scale_colour_sir()`](https://amr-for-r.org/reference/plot.html).

- upset_grid:

  Logical indicating whether to show marker combinations as an upset
  plot-style grid (default `FALSE`, so that each row is instead labelled
  with a printed list of markers).

- marker_label_space:

  Relative width of plotting area to provide to the marker list/grid.
  (Default `NULL`, which results in a default value of 3 when
  `upset_grid=FALSE` and 1 otherwise).

- plot_category:

  Logical indicating whether to include a stacked bar plot showing, for
  each marker combination, the proportion of samples with each phenotype
  classification (specified by the `pheno` column in the input file).
  Default is `TRUE`.

- print_category_counts:

  Logical indicating whether, if `plot_category=TRUE`, to print the
  number of strains in each resistance category for each marker
  combination in the plot. Default is `FALSE`.

- plot_ppv:

  Logical indicating whether to plot the estimates for positive
  predictive value, for each marker combination (default `TRUE`).

- plot_assay:

  Logical indicating whether to plot the distribution of MIC/disk assay
  values, for each marker combination (default `FALSE`).

- assay:

  A character string indicating whether to plot MIC or disk diffusion
  data. Must be one of:

  - `NULL`: (default) if no assay data is to be plotted

  - `"mic"`: plot MIC data stored in column `mic`

  - `"disk"`: plot disk diffusion data stored in column `disk`

- boxplot_col:

  Colour for lines of the box plots summarising the MIC distribution for
  each marker combination. Default is `"grey"`.

- species:

  Optional. Species name used to retrieve clinical breakpoints for
  annotation of the assay distribution plot.

- bp_site:

  Optional. Breakpoint site (e.g. "Non-meningitis") used when retrieving
  clinical breakpoints.

- guideline:

  Guideline used for breakpoint lookup. Default is `"EUCAST 2025"`.

- bp_S:

  (optional) S breakpoint to add to assay distribution plot (numerical).

- bp_R:

  (optional) R breakpoint to add to assay distribution plot (numerical).

- ecoff_bp:

  (optional) ECOFF breakpoint to add to assay distribution plot
  (numerical).

- pd:

  A
  [`ggplot2::position_dodge()`](https://ggplot2.tidyverse.org/reference/position_dodge.html)
  object controlling horizontal spacing of points and confidence
  intervals in the PPV plot. Default is `position_dodge(width = 0.8)`.

## Value

A list containing the following elements:

- `plot`: A grid of the requested plots

- `binary_matrix`: A copy of the genotype-phenotype binary matrix
  (either provided as input or generated by the function)

- `summary`: A data frame summarizing each marker combination observed,
  including number of resistant isolates, positive predictive values,
  and median assay values (and interquartile range) where relevant.

## Examples

``` r
if (FALSE) { # \dontrun{
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

# Generate binary matrix
binary_matrix <- get_binary_matrix(
  geno_table = ecoli_geno,
  pheno_table = ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "pheno_clsi",
  keep_assay_values = TRUE,
  keep_assay_values_from = "mic"
)

# Run ppv analysis using this binary_matrix
ppv(binary_matrix)

# Alternatively, generate binary matrix and run ppv() in one step
ppv(
  geno_table = ecoli_geno,
  pheno_table = ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "pheno_clsi"
)
} # }
```
