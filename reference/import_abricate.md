# Import and Process ABRicate Results

This function imports and processes ABRicate results, extracting
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

  A character string specifying a dataframe or path to the ABRicate
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

  A character string specifying which AMR gene database ABRicate was run
  with (default `"resfinder"`; `"ncbi"` is also supported).

## Value

A data frame with the processed genotype data, with harmonised gene
names, mapped drug agents, and drug classes which can be used for other
functions of the ARMgen package:

- `id`: The sample identifier (`character`).

- `marker`: The name of the genotype marker, as it appears in the `GENE`
  column of the input file (`character`).

- `gene`: The name of the gene product, as it appears in the `PRODUCT`
  column of the input file (`character`).

- `drug_class`: Name of the antibiotic group associated with the
  genotype marker, compatible with AMR pkg, parsed from the `RESISTANCE`
  column of the input file which depends on the database that ABRicate
  was run with (`character`).

- `drug`: Name of the specific antibiotic associated with the genotype
  marker, compatible with AMR pkg, parsed from the `RESISTANCE` column
  of the input file (`ab`). Value `NA` is assigned when the markers are
  annotated with a class only and not a specific antibiotic.

- `variation type`: Type of variation, i.e. `Gene presence detected`, as
  ABRicate only detects presence/absence of genes in the query database.
  ... Other fields specific to the input file

## Details

The function performs the following steps:

- Reads the ABRicate output table via the internal `process_input`
  function.

- Standardises the sample column name 'id'.

- Assigns standardised column names for genes, markers, and sets
  variation type to "Gene presence detected".

- Splits multiple resistance annotations (separated by semicolons) into
  separate rows.

- Converts drug agent names and classes to terms recognised by the AMR
  package.

## Examples

``` r
if (FALSE) { # \dontrun{
geno_table <- import_abricate("path/to/abricate_resfinder.tsv")

geno_table2 <- import_abricate("path/to/abricate_ncbi.tsv", db = "ncbi")
} # }
```
