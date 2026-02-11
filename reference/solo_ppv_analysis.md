# Perform Solo PPV Analysis for AMR Markers

This function performs a Positive Predictive Value (PPV) analysis for
AMR markers associated with a given antibiotic and drug class. It
calculates the PPV for solo markers and visualizes the results using
various plots.

## Usage

``` r
solo_ppv_analysis(
  geno_table,
  pheno_table,
  antibiotic = NULL,
  drug_class_list = NULL,
  geno_sample_col = NULL,
  pheno_sample_col = NULL,
  sir_col = NULL,
  ecoff_col = "ecoff",
  icat = FALSE,
  marker_col = "marker",
  binary_matrix = NULL,
  min = 1,
  axis_label_size = 9,
  pd = position_dodge(width = 0.8),
  excludeRanges = NULL,
  colours_SIR = c(S = "#3CAEA3", SDD = "#8FD6C4", I = "#F6D55C", R = "#ED553B"),
  colours_ppv = c(R = "maroon", I = "skyblue", NWT = "navy")
)
```

## Arguments

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
  specifying the antibiotic of interest to filter phenotype data. The
  value must match one of the entries in the `drug_agent` column of
  `pheno_table`. Only used if `binary_matrix` not provided or if
  breakpoints required.

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

- icat:

  A logical indicating whether to calculate PPV for `"I"` (if such a
  category exists in the phenotype column) (default `FALSE`).

- marker_col:

  A character string specifying the column name in `geno_table`
  containing the marker identifiers. Default `"marker"`. Only used if
  `binary_matrix` not provided.

- binary_matrix:

  A data frame containing the original binary matrix output from the
  [`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md)
  function. If not provided (or set to `NULL`), user must specify
  `geno_table`, `pheno_table`, `antibiotic`, `drug_class_list` and
  optionally `geno_sample_col`, `pheno_sample_col`, `sir_col`,
  `ecoff_col`, `marker_col` to pass to
  [`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md).

- min:

  Minimum number of genomes with the solo marker, to include the marker
  in the plot (default `1`).

- axis_label_size:

  Font size for axis labels in the PPV plot (default `9`).

- pd:

  A
  [`ggplot2::position_dodge()`](https://ggplot2.tidyverse.org/reference/position_dodge.html)
  object controlling horizontal spacing of points and confidence
  intervals in the PPV plot. Default `position_dodge(width = 0.8)`.

- excludeRanges:

  Vector of phenotype categories (comprised of `"R"`, `"I"`, `"NWT"`)
  for which we should ignore MIC values expressed as ranges when
  calculating PPVs. To include MICs expressed as ranges set this to
  `NULL`.

- colours_SIR:

  A named vector of colours for the percentage bar plot. The names
  should be the phenotype categories (e.g., `"R"`, `"I"`, `"S"`), and
  the values should be valid colour names or hexadecimal colour codes.
  Default values are those used in the AMR package
  [`AMR::scale_fill_sir()`](https://amr-for-r.org/reference/plot.html).

- colours_ppv:

  A named vector of colours for the plot of PPV estimates. The names
  should be `"R"`, `"I"`, `"NWT"`, and the values should be valid colour
  names or hexadecimal colour codes.

## Value

A list containing the following elements:

- `solo_stats`: A data frame summarizing the PPV for resistance (R vs
  S/I) and NWT (R/I vs S), including the number of positive hits, sample
  size, PPV, and 95% confidence intervals for each marker.

- `combined_plot`: A combined ggplot object showing the PPV plot for the
  solo markers, and a bar plot for the phenotype distribution.

- `solo_binary`: A dataframe with binary values indicating the presence
  or absence of the solo markers.

- `solo_binary_norange`: A dataframe with binary values indicating the
  presence or absence of the solo markers, excluding samples for which
  MIC or disk measures are expressed as ranges.

- `amr_binary`: A copy of the genotype-phenotype binary matrix for all
  markers (either provided as input or generated by the function)

- `plot_order`: Ordered list of rows in the output plot (to facilitate
  alignment with other plots)

## Details

The function analyzes the predictive power of individual AMR markers
when they are found 'solo' in the genome with no other markers
associated with the same class. The phenotype data are matched with
genotype presence/absence and then stratified to compute PPV for
resistance and non-wild-type interpretations. The function also
generates plots to aid in interpretation.

## Examples

``` r
if (FALSE) { # \dontrun{
geno_table <- import_amrfp(ecoli_geno_raw, "Name")
head(ecoli_ast)

# Generate binary matrix
binary_matrix <- get_binary_matrix(
  geno_table = geno_table,
  pheno_table = ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "pheno_clsi",
  keep_assay_values = TRUE,
  keep_assay_values_from = "mic"
)

# Run solo PPV analysis plot analysis using this binary_matrix
solo_ppv_analysis(binary_matrix = binary_matrix)

# Alternatively, generate binary matrix and solo PPV analysis in one step
# (note `antibiotic` and `drug_class` list are optional here, and only used
# for titling the plot)
soloPPV_cipro <- solo_ppv_analysis(
  binary_matrix = binary_matrix,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones")
)

soloPPV_cipro$solo_stats
soloPPV_cipro$combined_plot
} # }
```
