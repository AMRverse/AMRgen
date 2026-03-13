# Import and Process Abricate Results

This function imports and processes Abricate results, extracting
antimicrobial resistance (AMR) elements and mapping them to standardised
antibiotic names and drug classes. Currently supports results generated
using the ResFinder database.

## Usage

``` r
import_abricate(
  input_table,
  sample_col = "FILE",
  gene_col = "GENE",
  product_col = "PRODUCT",
  ab_col = "RESISTANCE",
  db = "resfinder"
)
```

## Arguments

- input_table:

  A character string specifying a dataframe or path to the Abricate
  results table.

- sample_col:

  A character string specifying the column that identifies samples in
  the dataset (default `"FILE"`).

- gene_col:

  A character string specifying the column that identifies gene symbols
  in the dataset (default `"GENE"`).

- product_col:

  A character string specifying the column that identifies product names
  in the dataset (default `"PRODUCT"`).

- ab_col:

  A character string specifying the column that identifies which drug/s
  each detected gene is associated with (default `"RESISTANCE"`).

- db:

  A character string specifying which AMR gene database Abricate was run
  with (currently only `"resfinder"` is supported).

## Value

A tibble containing the processed AMR elements, with harmonised gene
names, mapped drug agents, and drug classes which can be used for other
functions of the ARMgen package.

## Details

The function performs the following steps:

- Reads the Abricate output table via the internal `process_input`
  function.

- Standardises the sample column name 'id'.

- Assigns standardised column names for genes, markers, and sets
  variation type to "Gene presence detected".

- Splits multiple resistance annotations (separated by semicolons) into
  separate rows.

- Converts drug agent names to the `"ab"` class from the AMR package and
  maps these to classes compatible with the output of
  [`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
  and
  [`import_kleborate()`](https://AMRverse.github.io/AMRgen/reference/import_kleborate.md).

## Examples

``` r
if (FALSE) { # \dontrun{
geno_table <- import_abricate("path/to/abricate_resfinder.tsv")
} # }
```
