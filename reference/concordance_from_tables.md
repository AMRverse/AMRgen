# Assess concordance between tables of observed and predicted phenotypes

Calculates concordance from two data frames in standard phenotype table
format, one with observed phenotype data and the other with predicted
phenotype calls (e.g. output from AMRrules).

## Usage

``` r
concordance_from_tables(
  pheno_table,
  pheno_pred_table,
  true_SIR_col = "pheno_eucast",
  true_ecoff_col = "ecoff",
  pred_SIR = "clinical category",
  pred_ecoff = "phenotype",
  sample_col = "id",
  drug_col = "drug_agent",
  measure_col = "mic"
)
```

## Arguments

- pheno_table:

  A tibble or data frame containing real (observed) phenotype data, in
  the format output by
  [`import_pheno()`](https://AMRverse.github.io/AMRgen/reference/import_pheno.md).

- pheno_pred_table:

  A tibble or data frame containing predicted phenotype calls, e.g.
  AMRrules predictions imported in AMRgen using
  [`import_amrrules_predictions()`](https://AMRverse.github.io/AMRgen/reference/import_amrrules_predictions.md).

- true_SIR_col:

  Character. Name of the column containing S/I/R calls interpreted from
  real data against clinical breakpoints (default `"pheno_eucast"`).

- true_ecoff_col:

  Character. Name of the column containing WT/NWT calls interpreted from
  real data against ECOFF (default `"ecoff"`).

- pred_SIR:

  Character. Name of the column containing S/I/R calls predicted from
  genotypes (default `"clinical category"`).

- pred_ecoff:

  Character. Name of the column containing WT/NWT calls predicted from
  genotypes (default `"phenotype"`).

- sample_col:

  Character. Name of the column containing sample identifiers. This must
  be the same in both tables (default `"id"`).

- drug_col:

  Character. Name of the column containing drug agent identifiers. This
  must be the same in both tables (default `"drug_agent"`).

- measure_col:

  Character. Name of the column containing observed MIC or disk
  measurements for plotting. Valid options are `"mic"` or `"disk"`
  (default `"mic"`).

## Value

A named list with the following elements:

- obs_pred:

  A copy of the pheno_table with the predictions from pheno_pred_table
  merged in.

- metrics:

  A tibble listing concordance metrics, comparing observed vs predicted
  R and/or NWT calls, for all drugs with matched samples in both input
  tables.

- drugs_R:

  A long form tibble listing the number of samples with each combination
  of observed and predicted binary calls for R, for all drugs with
  matched samples in both input tables.

- drugs_NWT:

  A long form tibble listing the number of samples with each combination
  of observed and predicted binary calls for NWT, for all drugs with
  matched samples in both input tables.

- plot_R:

  A list containing two-panel plots summarising observed vs predicted
  clinical categories, one for each drug.

- plot_NWT:

  A list containing two-panel plots summarising observed vs predicted
  WT/NWT categories, one for each drug.

- plot_R_dist:

  A list containing plots of the assay value distributions (MIC or disk
  zones), coloured by S/I/R prediction, one for each drug.

- plot_NWT_dist:

  A list containing plots of the assay value distributions (MIC or disk
  zones), coloured by WT/NWT prediction, one for each drug.

## Examples

``` r
if (FALSE) { # \dontrun{
concordance_from_tables(true_phenotypes, predicted_phenotypes)
} # }
```
