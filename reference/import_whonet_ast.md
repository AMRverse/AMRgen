# Import and Process AST Data from WHONET Output Files

This function imports AST data from WHONET software output files (wide
CSV format) and converts it to the standardised long-format used by
AMRgen.

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
