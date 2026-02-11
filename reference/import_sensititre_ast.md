# Import and Process AST Data from Sensititre Output Files

This function imports AST data from Sensititre instrument output files
(UTF-16LE encoded, tab-separated, no header row) and converts it to the
standardised long-format used by AMRgen.

## Usage

``` r
import_sensititre_ast(
  input,
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

  Path to a Sensititre output text file

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
