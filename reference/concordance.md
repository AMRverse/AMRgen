# Calculate genotype-phenotype concordance from binary matrix

Compares genotypes (presence of resistance markers) to observed
phenotypes (R vs S, and/or NWT vs NWT) using a binary matrix from
[`get_binary_matrix()`](https://amrgen.org/reference/get_binary_matrix.md).
A genotypic prediction variable is defined on the basis of genotype
marker data (based on a variety of possible rules including
any/all/minimum number of markers; specifically those markers or
combinations exceeding a threshold positive predictive value (PPV);
predictions from a logistic regression model; or a user-defined field
providing predictions). This genotypic prediction is then compared to
the observed phenotypes using standard classification metrics (via the
`yardstick` pkg) and AMR-specific error rates (major error, ME and very
major error, VME) per ISO 20776-2 (and see [FDA
definitions](https://www.fda.gov/medical-devices/guidance-documents-medical-devices-and-radiation-emitting-products/antimicrobial-susceptibility-test-ast-systems-class-ii-special-controls-guidance-industry-and-fda)).
Supports evaluating both R and NWT outcomes in a single call, with
flexible prediction rules and marker inclusion options.

## Usage

``` r
concordance(
  binary_matrix,
  markers = NULL,
  exclude_markers = NULL,
  ppv_threshold = NULL,
  solo_ppv_results = NULL,
  ppv_results = NULL,
  prediction_col = NULL,
  truth = c("R", "NWT"),
  prediction_rule = "any",
  min_count = NULL,
  logreg_results = NULL,
  pval_threshold = NULL
)
```

## Arguments

- binary_matrix:

  A data frame output by
  [`get_binary_matrix()`](https://amrgen.org/reference/get_binary_matrix.md),
  containing one row per sample, columns indicating binary phenotypes
  (`R`, `I`, `NWT`) and binary marker presence/absence.

- markers:

  A character vector of marker column names to include in a new binary
  outcome variable `genotypic_prediction`. Default `NULL` includes all
  marker columns.

- exclude_markers:

  A character vector of marker column names to exclude from the
  genotypic prediction. Applied after `markers` filtering.

- ppv_threshold:

  A numeric PPV threshold (0-1). Used for solo PPV-based marker
  filtering when `solo_ppv_results` is provided, or as the combination
  PPV threshold when `prediction_rule = "combo_ppv"`.

- solo_ppv_results:

  Output of
  [`solo_ppv_analysis()`](https://amrgen.org/reference/solo_ppv_analysis.md),
  used for PPV-based marker filtering when `ppv_threshold` is set.

- ppv_results:

  Output of [`ppv()`](https://amrgen.org/reference/ppv.md), required
  when `prediction_rule = "combo_ppv"`. The `summary` table from this
  object is used to identify marker combinations with PPV \>=
  `ppv_threshold` for the relevant outcome (`R.ppv` or `NWT.ppv`).
  Samples whose marker combination matches any passing combination are
  predicted positive.

- prediction_col:

  A character string naming a column in `binary_matrix` that contains a
  user-defined prediction (coded 0/1). When supplied, all marker
  filtering and prediction generation are bypassed; the specified column
  is used directly as the prediction for all outcomes in `truth`. This
  allows arbitrary prediction logic to be evaluated using the
  concordance metrics.

- truth:

  A character vector specifying the phenotypic truth column(s) to
  evaluate: `"R"` (resistant vs susceptible/intermediate), `"NWT"`
  (non-wildtype vs wildtype), or `c("R", "NWT")` (default) to evaluate
  both.

- prediction_rule:

  The rule for generating genotypic predictions: `"any"` (default)
  predicts positive if any marker is present (after applying filters as
  specified by `markers`, `exclude_markers`, `ppv_threshold`,
  `pval_threshold`, `min_count`); `"all"` predicts positive only if all
  markers are present (after applying filters as specified by `markers`,
  `exclude_markers`, `ppv_threshold`, `pval_threshold`, `min_count`); a
  positive integer predicts positive if at least that many markers are
  present (after applying filters as specified by `markers`,
  `exclude_markers`, `ppv_threshold`, `pval_threshold`, `min_count`);
  `"logistic"` uses a logistic regression model from `logreg_results` to
  predict outcomes for each sample; `"combo_ppv"` predicts positive if a
  sample's marker combination has PPV \>= `ppv_threshold` in
  `ppv_results`.

- min_count:

  An integer or `NULL`. Exclude markers with total frequency (column sum
  in binary_matrix) below this value. Default `NULL` (no filtering).

- logreg_results:

  Output of
  [`amr_logistic()`](https://amrgen.org/reference/amr_logistic.md). Used
  for p-value filtering (when `pval_threshold` is set) and for
  `prediction_rule = "logistic"`.

- pval_threshold:

  A numeric p-value threshold. Exclude markers with logistic regression
  p-value \>= this value. Requires `logreg_results`.

## Value

An S3 object of class `"amr_concordance"`, a list containing:

- `conf_mat`: Named list of yardstick confusion matrix objects (e.g.
  `list(R = <cm>, NWT = <cm>)`).

- `metrics`: A tibble with columns `outcome`, `metric`, `estimate`.

- `data`: The input binary matrix with added `R_pred` and/or `NWT_pred`
  columns.

- `markers_used`: Named list of character vectors of markers used per
  outcome.

- `truth_col`: Character vector of truth columns evaluated.

- `n`: Named integer vector of sample counts per outcome.

- `prediction_rule`: The prediction rule used.

## Details

The function identifies marker columns as all columns not in the
reserved set (`id`, `pheno`, `ecoff`, `R`, `I`, `NWT`, `mic`, `disk`).
It then applies filtering in order: inclusion list `markers`, exclusion
list `exclude_markers`, `min_count`, `ppv_threshold` filtering, then
`pval_threshold` filtering.

Marker filtering is performed per outcome (R, NWT), since PPV-based and
p-value-based filters are category/outcome-specific.

When `prediction_col` is supplied, all marker filtering and prediction
generation are bypassed. The named column (0/1-coded) is used directly
as the genotypic prediction. This is useful when the user has computed a
custom prediction (e.g. based on marker combinations from
[`ppv()`](https://amrgen.org/reference/ppv.md)) and wants to evaluate
concordance metrics against the truth columns.

When `prediction_rule = "combo_ppv"`, the function calls
[`get_combo_matrix()`](https://amrgen.org/reference/get_combo_matrix.md)
internally to derive a combination identifier for each sample. Samples
whose combination identifier matches any entry in `ppv_results$summary`
with outcome PPV \>= `ppv_threshold` are predicted positive. The unique
individual markers contributing to passing combinations are reported in
`markers_used`.

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

[`get_binary_matrix()`](https://amrgen.org/reference/get_binary_matrix.md),
[`solo_ppv_analysis()`](https://amrgen.org/reference/solo_ppv_analysis.md),
[`amr_logistic()`](https://amrgen.org/reference/amr_logistic.md),
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

# Basic concordance for both R and NWT
result <- concordance(binary_matrix)
result

# Single outcome
result <- concordance(binary_matrix, truth = "R")

# Exclude specific markers
result <- concordance(binary_matrix, exclude_markers = c("qnrS1"))

# Filter markers by solo PPV threshold
solo_ppv <- solo_ppv_analysis(binary_matrix = binary_matrix)
result <- concordance(
  binary_matrix,
  ppv_threshold = 0.5,
  solo_ppv_results = solo_ppv
)

# Require at least 2 markers present for prediction
result <- concordance(binary_matrix, prediction_rule = 2)

# Use logistic regression model for prediction
logreg <- amr_logistic(binary_matrix = binary_matrix)
result <- concordance(
  binary_matrix,
  prediction_rule = "logistic",
  logreg_results = logreg
)

# Predict based on marker combinations with PPV >= 0.5 (from ppv())
ppv_res <- ppv(binary_matrix = binary_matrix)
result <- concordance(
  binary_matrix,
  prediction_rule = "combo_ppv",
  ppv_results = ppv_res,
  ppv_threshold = 0.5
)

# Use a custom user-defined prediction column
binary_matrix$my_pred <- as.integer(binary_matrix$gyrA_S83L == 1 | binary_matrix$gyrA_D87N == 1)
result <- concordance(binary_matrix, prediction_col = "my_pred", truth = "R")

# Access components
result$conf_mat$R
result$metrics
result$markers_used

# predictions vs observed SIR calls
result$data %>% count(R_pred, pheno)
} # }
```
