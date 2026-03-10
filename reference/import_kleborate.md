# Import and Process Kleborate Results

This function imports and processes genotyping results from Kleborate
(https://github.com/klebgenomics/Kleborate), extracting antimicrobial
resistance determinants and mapping them to standardised drug classes.

## Usage

``` r
import_kleborate(input_table, sample_col = "strain")
```

## Arguments

- input_table:

  A character string specifying a dataframe or path to the Kleborate
  results table (TSV format).

- sample_col:

  A character string specifying the column that identifies samples in
  the dataset (default `strain`).

## Value

A tibble containing the processed AMR determinants and drug classes that
is AMRgen compatible.

## Details

The function performs the following steps:

- Reads the Kleborate output table.

- Transforms Kleborate output into long form (i.e., one AMR determinant
  per row).

- Maps Kleborate drug classes to standardised drug class names. This
  processing ensures compatibility with downstream AMRgen analysis
  workflows.

## Examples

``` r
if (FALSE) { # \dontrun{
# example Kleborate data from EUSCAPE project
kleborate_raw

# import first few rows of this data frame and parse it as AMRfp data
kleborate_geno <- import_kleborate(kleborate_raw %>% head(n = 10), "strain")
geno
} # }
```
