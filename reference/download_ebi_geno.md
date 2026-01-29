# Download antimicrobial genotype data from the EBI AMR Portal

This function will retrieve genotype data from the EBI AMR Portal,
https://www.ebi.ac.uk/amr. The portal uses AMRfinderplus to identify
AMR-associated genotypes, but the results are processed and not all
fields returned by AMRfinderplus are included. See
https://www.ebi.ac.uk/amr/about/#AMR-Genotypes for more information, and
https://github.com/ncbi/amr/wiki/class-subclass for valid class and
subclass terms.

## Usage

``` r
download_ebi_geno(
  user_genus = NULL,
  user_subclass = NULL,
  user_class = NULL,
  user_antibiotic_name = NULL,
  user_release = NULL
)
```

## Arguments

- user_genus:

  String specifying a bacterial genus to return data for (default NULL,
  will pull all taxa)

- user_subclass:

  String specifying an antibiotic subclass to filter on (default NULL,
  check NCBI Subclass for valid terms).

- user_class:

  String specifying an antibiotic subclass to filter on (default NULL,
  check NCBI Class for valid terms).

- user_antibiotic_name:

  String specifying an antibiotic to return data for (default NULL, will
  pull all antibiotics).

- user_release:

  String specifying the data release to download (default NULL, will
  pull latest release).

## Value

A data frame containing EBI genotype data

## Examples

``` r
if (FALSE) { # \dontrun{
ebi_genotypes <- import_ebi()

amp_sal_ebi <- import_ebi(
    user_genus="Salmonella",
    user_release="2025-12",
    user_subclass="beta-lactam"
)
} # }
```
