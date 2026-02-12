# Import and process antimicrobial phenotype data exported from MicroScan instruments

This function imports antimicrobial susceptibility testing (AST) data
from MicroScan instrument output files (wide CSV format) and converts it
to the standardised long-format used by AMRgen. Supports English,
Spanish, French, German, and Portuguese column names (auto-detected from
metadata columns).

## Usage

``` r
import_microscan_ast(
  input,
  sample_col = NULL,
  source = NULL,
  species = NULL,
  ab = NULL,
  instrument_guideline = NULL,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE
)
```

## Arguments

- input:

  A dataframe or path to a CSV/TSV file containing MicroScan AST output
  data

- sample_col:

  Column name for sample identifiers. Default: NULL (auto-detected from
  language-specific column names)

- source:

  Optional source value to record for all data points

- species:

  Optional species override for phenotype interpretation

- ab:

  Optional antibiotic override for phenotype interpretation

- instrument_guideline:

  Optional guideline used by the instrument for SIR interpretation

- interpret_eucast:

  Interpret against EUCAST breakpoints

- interpret_clsi:

  Interpret against CLSI breakpoints

- interpret_ecoff:

  Interpret against ECOFF values

## Value

Standardised AST data frame
