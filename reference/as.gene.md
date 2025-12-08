# Gene Class and AMR Parsing Functions

Functions to create a custom "gene" class and parse AMR data.

## Usage

``` r
as.gene(x)
```

## Arguments

- x:

  A character vector to be converted to a "gene" class.

## Value

For `as.gene`, an object of class `"gene"`.

## Examples

``` r
# Create a gene object
gene <- as.gene(c("gene1", "gene2"))
gene
#> An object of class 'gene':
#> [1] "gene1" "gene2"

if (FALSE) { # \dontrun{
# Parse AMR data
parsed_data <- import_amrfp("path/to/input_table.tsv", "SampleID")
} # }
```
