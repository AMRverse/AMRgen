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
  antibiotic,
  drug_class_list,
  geno_sample_col = NULL,
  pheno_sample_col = NULL,
  sir_col = NULL,
  ecoff_col = "ecoff",
  icat = FALSE,
  marker_col = "marker",
  keep_assay_values = TRUE,
  min = 1,
  axis_label_size = 9,
  pd = position_dodge(width = 0.8),
  plot_cols = c(R = "maroon", I = "skyblue", NWT = "navy"),
  excludeRanges = c("NWT")
)
```

## Arguments

- geno_table:

  A data frame containing genotype data, including at least one column
  labeled `drug_class` for drug class information and one column for
  sample identifiers (specified via `geno_sample_col` otherwise it is
  assumed the first column contains identifiers).

- pheno_table:

  A data frame containing phenotype data, which must include a column
  `drug_agent` (with the antibiotic information) and a column with the
  resistance interpretation (S/I/R, colname specified via `sir_col`).

- antibiotic:

  A character string specifying the antibiotic of interest to filter
  phenotype data. The value must match one of the entries in the
  `drug_agent` column of `pheno_table`.

- drug_class_list:

  A character vector of drug classes to filter genotype data for markers
  related to the specified antibiotic. Markers in `geno_table` will be
  filtered based on whether their `drug_class` matches any value in this
  list.

- geno_sample_col:

  A character string (optional) specifying the column name in
  `geno_table` containing sample identifiers. Defaults to `NULL`, in
  which case it is assumed the first column contains identifiers.

- pheno_sample_col:

  A character string (optional) specifying the column name in
  `pheno_table` containing sample identifiers. Defaults to `NULL`, in
  which case it is assumed the first column contains identifiers.

- sir_col:

  A character string specifying the column name in `pheno_table` that
  contains the resistance interpretation (SIR) data. The values should
  be interpretable as "R" (resistant), "I" (intermediate), or "S"
  (susceptible).

- ecoff_col:

  A character string (optional) specifying the column name in
  `pheno_table` that contains resistance interpretations (SIR) made
  against the ECOFF rather than a clinical breakpoint. The values should
  be interpretable as "R" (resistant), "I" (intermediate), or "S"
  (susceptible).

- icat:

  A logical indicating whether to calculate PPV for "I" (if such a
  category exists in the phenotype column) (default FALSE).

- marker_col:

  A character string specifying the column name in `geno_table`
  containing the marker identifiers. Defaults to `"marker"`.

- keep_assay_values:

  A logical indicating whether to include columns with the raw phenotype
  assay data in the binary matrix. Assumes there are columns labelled
  "mic" and/or "disk"; these will be added to the output table if
  present. Defaults to `TRUE`.

- min:

  Minimum number of genomes with the solo marker, to include the marker
  in the plot (default 1).

- axis_label_size:

  Font size for axis labels in the PPV plot (default 9).

- pd:

  Position dodge, i.e. spacing for the R/NWT values to be positioned
  above/below the line in the PPV plot. Default 'position_dodge(width =
  0.8)'.

- plot_cols:

  A named vector of colors for the plot. The names should be the
  phenotype categories (e.g., "R", "I", "S", "NWT"), and the values
  should be valid color names or hexadecimal color codes. Default colors
  are provided for resistant ("R"), intermediate ("I"), susceptible
  ("S"), and non-wild-type ("NWT").

- excludeRanges:

  Vector of phenotype categories (comprised of "R", "I", "NWT") for
  which we should ignore MIC values expressed as ranges when calculating
  PPVs. Default c("NWT"), as calling against ECOFF with the AMR package
  currently does not interpret ranges correctly. To include MICs
  expressed as ranges set this to NULL.

## Value

A list containing the following elements:

- `solo_stats`: A dataframe summarizing the PPV for resistance (R vs
  S/I) and NWT (R/I vs S), including the number of positive hits, sample
  size, PPV, and 95% confidence intervals for each marker.

- `combined_plot`: A combined ggplot object showing the PPV plot for the
  solo markers, and a bar plot for the phenotype distribution.

- `solo_binary`: A dataframe with binary values indicating the presence
  or absence of the solo markers.

- `amr_binary`: A dataframe with binary values for the AMR markers,
  based on the input genotype and phenotype data.

## Details

The function analyzes the predictive power of individual AMR markers
that belong to a specified drug class and are uniquely associated with
one class. The phenotype data are matched with genotype presence/absence
and then stratified to compute PPV for resistance and non-wild-type
interpretations. It also generates plots to aid in interpretation.

## Examples

``` r
if (FALSE) { # \dontrun{
geno_table <- import_amrfp(ecoli_geno_raw, "Name")
head(ecoli_ast)
soloPPV_cipro <- solo_ppv_analysis(
  geno_table = geno_table,
  pheno_table = ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype"
)
soloPPV_cipro$solo_stats
soloPPV_cipro$combined_plot
} # }
```
