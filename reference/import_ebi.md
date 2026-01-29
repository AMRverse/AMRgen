# Download antimicrobial phenotype data from the EBI AMR Portal

Download antimicrobial phenotype data from the EBI AMR Portal

## Usage

``` r
import_ebi(user_genus = NULL, user_release = NULL, user_antibiotic_name = NULL)
```

## Arguments

- user_genus:

  String specifying a bacterial genus to download data for (default
  NULL, will pull all taxa)

- user_release:

  String specifying the data release to download (default NULL, will
  pull latest release)

- user_antibiotic_name:

  String specifying an antibiotic to download data for (default NULL,
  will pull all antibiotics)

## Value

A data frame containing EBI genotype data

## Examples

``` r
if (FALSE) { # \dontrun{
import_ebi(
    user_genus="Salmonella",
    user_release="2025-12",
    user_antibiotic_name="ampicillin"
)
} # }
```
