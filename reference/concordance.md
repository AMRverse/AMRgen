# Calculate genotype-phenotype concordance from binary matrix

Compares genotypes (presence of resistance markers) to observed
phenotypes (resistant vs susceptible) using a binary matrix from
[`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md).
A genotypic prediction variable is defined on the basis of presence of
genotype markers (either any marker in the input table, or those defined
by an input inclusion list or exclusion list). This genotypic prediction
is then compared to the observed phenotypes using standard
classification metrics (via the `yardstick` pkg) and AMR-specific error
rates (major error, ME and very major error, VME) per ISO 20776-2 (and
see [FDA
definitions](https://www.fda.gov/medical-devices/guidance-documents-medical-devices-and-radiation-emitting-products/antimicrobial-susceptibility-test-ast-systems-class-ii-special-controls-guidance-industry-and-fda).

## Usage

``` r
concordance(
  binary_matrix,
  markers = NULL,
  exclude_markers = NULL,
  ppv_threshold = NULL,
  solo_ppv_results = NULL,
  truth = "R",
  prediction_rule = "any"
)
```

## Arguments

- binary_matrix:

  A data frame output by
  [`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md),
  containing one row per sample, columns indicating binary phenotypes
  (`R`, `I`, `NWT`) and binary marker presence/absence.

- markers:

  A character vector of marker column names to include in a summary
  binary outcome variable 'markers present'. Default `NULL` includes all
  marker columns.

- exclude_markers:

  A character vector of marker column names to exclude from the
  genotypic prediction. Applied after `markers` filtering.

- ppv_threshold:

  A numeric PPV threshold (0-1). Markers with solo PPV below this value
  are excluded. Requires `solo_ppv_results`.

- solo_ppv_results:

  Output of
  [`solo_ppv_analysis()`](https://AMRverse.github.io/AMRgen/reference/solo_ppv_analysis.md),
  used for PPV-based marker filtering when `ppv_threshold` is set.

- truth:

  A character string specifying the phenotype column to use: `"R"` (R vs
  S/I classifications, default) or `"NWT"` (nonwildtype vs wildtype).

- prediction_rule:

  A character string specifying the rule for generating genotypic
  predictions. Currently only `"any"` is supported: a sample is
  predicted positive if any included marker is present (value of 1).

## Value

An S3 object of class `"amr_concordance"`, a list containing:

- `conf_mat`: A yardstick confusion matrix object.

- `metrics`: A tibble with columns `.metric`, `.estimator`, `.estimate`
  containing sensitivity, specificity, ppv, npv, accuracy, kap, f_meas,
  VME, and ME.

- `data`: The input binary matrix with an added `geno_prediction`
  column.

- `markers_used`: Character vector of markers included in the
  prediction.

- `truth_col`: Which truth column was used (`"R"` or `"NWT"`).

- `n`: Total number of samples with non-missing truth values.

## Details

The function identifies marker columns as all columns not in the
reserved set (`id`, `pheno`, `ecoff`, `R`, `I`, `NWT`, `mic`, `disk`).
It then applies filtering in order: inclusion list `markers`, exclusion
list `exclude_markers`, then `ppv_threshold` filtering.

Genotypic prediction is generated per sample: if
`prediction_rule = "any"`, then a sample is predicted positive (1) if
any included marker equals 1.

Standard metrics (sensitivity, specificity, PPV, NPV, accuracy, kappa,
F-measure) are calculated using pkg `yardstick`. AMR-specific error
rates are computed internally:

- **VME** (Very Major Error): FN / (TP + FN) = 1 - sensitivity.
  Proportion of truly resistant isolates not predicted as such from
  genotype.

- **ME** (Major Error): FP / (TN + FP) = 1 - specificity. Proportion of
  truly susceptible isolates incorrectly predicted resistant from
  genotype.

## See also

[`get_binary_matrix()`](https://AMRverse.github.io/AMRgen/reference/get_binary_matrix.md),
[`solo_ppv_analysis()`](https://AMRverse.github.io/AMRgen/reference/solo_ppv_analysis.md),
[yardstick::yardstick](https://yardstick.tidymodels.org/reference/yardstick-package.html)

## Examples

``` r
if (FALSE) { # \dontrun{
geno_table <- import_amrfp(ecoli_geno_raw, "Name")

binary_matrix <- get_binary_matrix(
  geno_table = geno_table,
  pheno_table = ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "pheno_clsi"
)

# Basic concordance using all markers
result <- concordance(binary_matrix)
result

# Exclude specific markers
result <- concordance(binary_matrix, exclude_markers = c("qnrS1"))

# Filter markers by solo PPV threshold
solo_ppv <- solo_ppv_analysis(binary_matrix = binary_matrix)
result <- concordance(
  binary_matrix,
  ppv_threshold = 0.5,
  solo_ppv_results = solo_ppv
)

# Access components
result$conf_mat
result$metrics
result$markers_used
} # }
```
