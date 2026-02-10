# Add marker combinations to a binary geno-pheno matrix

Given a geno-pheno binary marker matrix, output by `get_binary_matrix`,
this function identifies marker combination, and adds combination ids,
and the count of markers detected per sample, as new columns.

## Usage

``` r
get_combo_matrix(binary_matrix, assay = NULL)
```

## Arguments

- binary_matrix:

  A geno-pheno binary matrix, output by `get_binary_matrix`

- assay:

  (optional) Name of an assay column to filter on, so that the matrix
  returned only includes samples with assay data of this type available

## Examples

``` r
if (FALSE) { # \dontrun{
combo_matrix <- get_combo_matrix(binary_matrix)
} # }
```
