# Co-occuring PPV Analysis

This function performs a Positive Predictive Value (PPV) analysis for
single or multiple AMR markers associated with a given antibiotic and
drug class. It calculates the PPV for single and co-occurring markers
and visualizes the results using UpSet-style matrix plot.

## Usage

``` r
cooccuring_ppv_analysis(
  geno_table,
  pheno_table,
  antibiotic,
  drug_class_list,
  geno_sample_col = NULL,
  pheno_sample_col = NULL,
  sir_col = NULL,
  marker_col = "marker",
  keep_assay_values = TRUE,
  min = 1,
  axis_label_size = 9,
  pd = position_dodge(width = 0.8),
  plot_cols = c(R = "IndianRed", NWT = "navy")
)
```

## Arguments

- geno_table:

  A data frame containing genotype data.

- pheno_table:

  A data frame containing phenotype data.

- antibiotic:

  String. The specific antibiotic to analyze (e.g., "Ceftazidime").

- drug_class_list:

  Vector of strings. The drug classes to filter for (e.g.,
  c("Cephalosporins")).

- geno_sample_col:

  String. Column name for sample IDs in genotype table.

- pheno_sample_col:

  String. Column name for sample IDs in phenotype table.

- sir_col:

  String. Column name containing S/I/R data in the phenotype table.

- marker_col:

  String. Default "marker". Name of the column identifying genes.

- keep_assay_values:

  Logical. Default TRUE. Whether to keep MIC/Disk columns.

- min:

  Integer. Default 1. Minimum number of isolates required for a profile
  to be plotted.

- axis_label_size:

  Integer. Font size for axis text.

- pd:

  Position dodge object for plotting.

- plot_cols:

  Named vector of colors for "R" and "NWT".

## Value

A list containing:

- `stats`: Data frame of calculated PPV statistics per profile.

- `plot`: The combined ggplot/patchwork object.

- `data`: The filtered data frame used for the analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
geno_table <- import_amrfp(ecoli_geno_raw, "Name")
head(ecoli_ast)
cooccuringPPV_cipro <- cooccuring_ppv_analysis(
  geno_table = geno_table,
  pheno_table = ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "pheno_provided"
)
} # }
```
