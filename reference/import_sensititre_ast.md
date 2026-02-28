# Import and process antimicrobial phenotype data exported from Sensititre instruments

This function imports antimicrobial susceptibility testing (AST) data
from Sensititre instrument output files (tab- or comma-separated,
optionally UTF-16LE encoded, no header row) and converts it to the
standardised long-format used by AMRgen.

## Usage

``` r
import_sensititre_ast(
  input,
  source = NULL,
  species = NULL,
  ab = NULL,
  instrument_guideline = NULL,
  id_col = 7,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE
)
```

## Arguments

- input:

  Path to a Sensititre output text file (tab- or comma-separated)

- source:

  Optional source value to record for all data points

- species:

  Optional species override for phenotype interpretation

- ab:

  Optional antibiotic override for phenotype interpretation

- instrument_guideline:

  Optional guideline used by the instrument for SIR interpretation

- id_col:

  Integer. Column index (1-based) of the sample identifier. Default is
  7, which corresponds to the sample accession column in standard
  Sensititre exports. Adjust if your file uses a different column for
  the sample ID (e.g. set to 2 for the plate/batch identifier column).

- interpret_eucast:

  Interpret against EUCAST breakpoints

- interpret_clsi:

  Interpret against CLSI breakpoints

- interpret_ecoff:

  Interpret against ECOFF values

## Value

Standardised AST data frame
