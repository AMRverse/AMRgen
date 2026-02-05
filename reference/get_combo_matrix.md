# Generate matrix of marker combinations

Given a geno-phenot binary marker matrix, output by `get_binary_matrix`,
this function constructs identifies marker combination, reshapes the
data to long format, and computes frequencies of marker combinations and
individual markers.

## Usage

``` r
get_combo_matrix(binary_matrix, assay = "mic")
```

## Arguments

- binary_matrix:

  A geno-pheno binary matrix, output by `get_binary_matrix`

- assay:

  (optional) Name of an assay column to filter on, so that the matrix
  returned only includes samples with assay data of this type available

## Value

A named list with three elements:

- combination_matrix:

  A long-format data frame with one row per sample–marker combination,
  including the marker presence/absence and a combination identifier.

- combination_freq:

  A data frame summarising the frequency and percentage of each unique
  marker combination.

- marker_freq:

  A data frame giving the prevalence (count) of each individual marker
  across all samples.

## Examples

``` r
if (FALSE) { # \dontrun{
combo_matrix <- get_combo_matrix(binary_matrix)
} # }
```
