# Generate Upset Plot

This function generates an upset plot showing summaries of phenotype
results (assay distributions, phenotype category percentages, and/or
predictive value for phenotype) for each combination of markers observed
in the data.

## Usage

``` r
amr_upset(
  binary_matrix = NULL,
  assay = "mic",
  min_set_size = 2,
  order = "value",
  geno_table,
  pheno_table,
  antibiotic = NULL,
  drug_class_list,
  geno_sample_col = NULL,
  pheno_sample_col = NULL,
  sir_col = NULL,
  ecoff_col = "ecoff",
  marker_col = "marker",
  plot_set_size = FALSE,
  plot_category = TRUE,
  print_category_counts = FALSE,
  print_set_size = FALSE,
  boxplot_col = "grey",
  SIR_col = c(S = "#3CAEA3", I = "#F6D55C", R = "#ED553B"),
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
  [`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md)
  function. If not provided (or set to `NULL`), user must specify
  `geno_table`, `pheno_table`, `antibiotic`, `drug_class_list` and
  optionally `geno_sample_col`, `pheno_sample_col`, `sir_col`,
  `ecoff_col`, `marker_col` to pass to
  [`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md).

- assay:

  A character string indicating whether to plot MIC or disk diffusion
  data. Must be one of:

  - `"mic"`: plot MIC data stored in column `mic`

  - `"disk"`: plot disk diffusion data stored in column `disk`

- min_set_size:

  An integer specifying the minimum size for a gene set to be included
  in the analysis and plots. Default is 2. Only marker combinations with
  at least this number of occurrences are included in the plots.

- order:

  A character string indicating the order of the combinations on the
  x-axis. Options are:

  - `""`: decreasing frequency of combinations

  - `"genes"`: order by the number of genes in each combination

  - `"value" (default)`: order by the median assay value (MIC or disk
    zone) for each combination

  - `"ppv"`: order by the PPV estimated for each combination

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

  Optional. Antibiotic name used to retrieve clinical breakpoints for
  annotation of the assay distribution plot.

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

- plot_set_size:

  Logical indicating whether to include a bar plot showing the set size
  (i.e., number of times each combination of markers is observed).
  Default is `FALSE`.

- plot_category:

  Logical indicating whether to include a stacked bar plot showing, for
  each marker combination, the proportion of samples with each phenotype
  classification (specified by the `pheno` column in the input file).
  Default is `TRUE`.

- print_category_counts:

  Logical indicating whether, if `plot_category=TRUE`, to print the
  number of strains in each resistance category for each marker
  combination in the plot. Default is `FALSE`.

- print_set_size:

  Logical indicating whether, if `plot_set_size=TRUE`, to print the
  number of strains with each marker combination on the plot. Default is
  `FALSE`.

- boxplot_col:

  Colour for lines of the box plots summarising the MIC distribution for
  each marker combination. Default is `"grey"`.

- SIR_col:

  A named vector of colours for the percentage bar plot. The names
  should be the phenotype categories (e.g., `"R"`, `"I"`, `"S"`), and
  the values should be valid color names or hexadecimal color codes.
  Default values are those used in the AMR package
  [`AMR::scale_colour_sir()`](https://amr-for-r.org/reference/plot.html).

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

## Value

A list containing the following elements:

- `plot`: A grid of plots displaying: (i) grid showing the marker
  combinations observed, MIC distribution per marker combination,
  frequency per marker and (optionally) phenotype classification and/or
  number of samples for each marker combination.

- `binary_matrix`: A copy of the genotype-phenotype binary matrix
  (either provided as input or generated by the function)

- `summary`: A data frame summarizing each marker combination observed,
  including median MIC (and interquartile range), number of resistant
  isolates, and positive predictive value for resistance.

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

# Run upset plot analysis using this binary_matrix
amr_upset(binary_matrix, assay = "mic")

# Alternatively, generate binary matrix and run ppv() in one step
amr_upset(
  assay = "mic",
  geno_table = ecoli_geno,
  pheno_table = ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "pheno_clsi"
)
} # }
```
