# Download antimicrobial genotype data from the EBI AMR Portal

This function will retrieve genotype data from the EBI AMR Portal,
https://www.ebi.ac.uk/amr. The portal uses AMRfinderplus to identify
AMR-associated genotypes, but the results are processed and not all
fields returned by AMRfinderplus are included. See
https://www.ebi.ac.uk/amr/about/#AMR-Genotypes for more information.

## Usage

``` r
download_ebi_geno(
  user_genus = NULL,
  user_release = NULL,
  user_antibiotic_name = NULL
)
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
  will pull all antibiotics).

## Value

A data frame containing EBI genotype data

## Examples

``` r
if (FALSE) { # \dontrun{
amp_sal_ebi <- import_ebi(
    user_genus="Salmonella",
    user_release="2025-12",
    user_antibiotic_name="ampicillin"
)
} # }
```
