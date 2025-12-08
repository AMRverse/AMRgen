# hamronize_data

hamronize_data

## Usage

``` r
harmonize_data(
  user_software_name,
  user_software_version,
  user_database_version,
  user_input_filename
)
```

## Arguments

- user_software_name:

  the analysis software used to screen genome data for AMR determinants
  (must be amrfinderplus, rgi, or resfinder)

- user_software_version:

  the version of the analysis software used to screen genome data

- user_database_version:

  the version of the database used

- user_input_filename:

  the name of the genotypic AMR data file

## Value

A data frame containing 'harmonized' AMR genotype data

## Examples

``` r
if (FALSE) { # \dontrun{
harmonize_data(
  user_software_name = "amrfinderplus",
  user_software_version = "3.12.8",
  user_input_filename = "ATB_Achromobacter_AFP.tsv",
  user_database_version = "2024-01-31.1"
)
} # }
```
