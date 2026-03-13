# Import and process antimicrobial phenotype data from WHONET files

This function imports antimicrobial susceptibility testing (AST) data
from WHONET software output files (wide CSV format) and converts it to
the standardised long-format used by AMRgen.

## Usage

``` r
import_whonet_ast(
  input,
  sample_col = "Identification number",
  source = NULL,
  species = NULL,
  ab = NULL,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE,
  include_patient_info = FALSE
)
```

## Arguments

- input:

  A dataframe or path to a CSV file containing WHONET AST output data

- sample_col:

  Column name for sample identifiers. Default: "Identification number"

- source:

  Optional source value to record for all data points

- species:

  Optional species override for phenotype interpretation

- ab:

  Optional antibiotic override for phenotype interpretation

- interpret_eucast:

  Interpret against EUCAST breakpoints

- interpret_clsi:

  Interpret against CLSI breakpoints

- interpret_ecoff:

  Interpret against ECOFF values

- include_patient_info:

  Include patient demographic columns in output

## Value

Standardised AST data frame

## Examples

``` r
# Built-in AMR package WHONET example dataset (standard uppercase format)
result <- import_whonet_ast(AMR::WHONET)
#> Warning: There were 2 warnings in `mutate()`.
#> The first warning was:
#> ℹ In argument: `mic = as.mic(...)`.
#> Caused by warning:
#> ! in `as.mic()`: 1311 results in column 'mic' truncated (9%) that were
#> invalid MICs: "S", "I", and "R"
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.
head(result)
#> # A tibble: 6 × 31
#>   id         drug_agent   mic  disk guideline method       platform disk_potency
#>   <chr>      <ab>       <mic> <dsk> <chr>     <chr>        <chr>    <chr>       
#> 1 fe41d7bafa AMP           NA    NA CLSI      disk diffus… NA       10          
#> 2 fe41d7bafa AMC           NA    NA EUCAST    disk diffus… NA       20          
#> 3 fe41d7bafa TZP           NA    NA EUCAST    disk diffus… NA       30          
#> 4 fe41d7bafa FEP           NA    NA EUCAST    disk diffus… NA       30          
#> 5 fe41d7bafa CTX           NA    NA EUCAST    disk diffus… NA       5           
#> 6 fe41d7bafa FOX           NA    NA EUCAST    disk diffus… NA       30          
#> # ℹ 23 more variables: pheno_provided <sir>, spp_pheno <mo>,
#> #   `Specimen number` <int>, Organism <chr>, Country <chr>, Laboratory <chr>,
#> #   collection_date <date>, `Specimen type` <chr>,
#> #   `Specimen type (Numeric)` <dbl>, Reason <chr>, `Isolate number` <int>,
#> #   `Organism type` <chr>, Serotype <chr>, `Beta-lactamase` <lgl>, ESBL <lgl>,
#> #   Carbapenemase <lgl>, `MRSA screening test` <lgl>,
#> #   `Inducible clindamycin resistance` <lgl>, Comment <chr>, …
if (FALSE) { # \dontrun{
# WHONET file with lowercase columns and _SIR suffix (e.g. amc_nd20, amc_nd20_SIR)
# import_whonet_ast("path/to/whonet_export.csv", sample_col = "patient_id")

# Include patient demographics in output
import_whonet_ast("path/to/whonet_export.csv",
  sample_col = "patient_id",
  include_patient_info = TRUE
)

# Interpret against EUCAST breakpoints for a specific species
result_interp <- import_whonet_ast(AMR::WHONET,
  interpret_eucast = TRUE,
  species = "Escherichia coli"
)
} # }
```
