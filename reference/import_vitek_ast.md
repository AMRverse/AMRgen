# Import and Process AST Data from Vitek Output Files

This function imports AST data from Vitek instrument output files (wide
CSV format) and converts it to the standardised long-format used by
AMRgen.

## Usage

``` r
import_vitek_ast(
  input,
  sample_col = "Lab ID",
  source = NULL,
  species = NULL,
  ab = NULL,
  instrument_guideline = NULL,
  use_expertized = TRUE,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE,
  include_dates = TRUE
)
```

## Arguments

- input:

  A dataframe or path to a CSV file containing Vitek AST output data

- sample_col:

  Column name for sample identifiers. Default: "Lab ID"

- source:

  Optional source value to record for all data points (e.g., dataset
  name or study identifier)

- species:

  Optional species override for phenotype interpretation

- ab:

  Optional antibiotic override for phenotype interpretation

- instrument_guideline:

  Optional guideline used by the Vitek instrument for SIR interpretation
  (e.g., "EUCAST 2025", "CLSI 2025"). Default: NULL

- use_expertized:

  Use expertized SIR (TRUE, default) or instrument SIR (FALSE)

- interpret_eucast:

  Interpret against EUCAST breakpoints

- interpret_clsi:

  Interpret against CLSI breakpoints

- interpret_ecoff:

  Interpret against ECOFF values

- include_dates:

  Include collection_date and testing_date in output

## Value

Standardised AST data frame
