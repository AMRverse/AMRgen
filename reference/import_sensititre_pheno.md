# Import and process antimicrobial phenotype data exported from Sensititre instruments

This function imports antimicrobial susceptibility testing (AST) data
from Sensititre instrument output files (tab- or comma-separated,
optionally UTF-16LE encoded, no header row) and converts it to the
standardised long-format used by AMRgen.

## Usage

``` r
import_sensititre_pheno(
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

  Optional string value to record in the `source` column for all data
  points (e.g., dataset name or study identifier)

- species:

  Optional string indicating a single species to use for phenotype
  interpretation (otherwise this is inferred per-sample from the input)

- ab:

  Optional string indicating a single antibiotic to use for phenotype
  interpretation (otherwise this is inferred per-sample from the input)

- instrument_guideline:

  Optional string indicating the guideline used by the instrument for
  SIR interpretation (e.g., "EUCAST 2025", "CLSI 2025"), used to record
  the `guideline` in the output file. Default: `NULL`

- id_col:

  Integer. Column index (1-based) of the sample identifier. Default is
  7, which corresponds to the sample accession column in standard
  Sensititre exports. Adjust if your file uses a different column for
  the sample ID (e.g. set to 2 for the plate/batch identifier column).

- interpret_eucast:

  A logical value (default is `FALSE`). If `TRUE`, the function will
  re-interpret the susceptibility phenotype (SIR) for each observation
  based on the MIC values, against EUCAST human breakpoints. These will
  be reported in a new column `pheno_eucast`, of class `sir`.

- interpret_clsi:

  A logical value (default is `FALSE`). If `TRUE`, the function will
  re-interpret the susceptibility phenotype (SIR) for each observation
  based on the MIC values, against CLSI human breakpoints. These will be
  reported in a new column `pheno_clsi`, of class `sir`.

- interpret_ecoff:

  A logical value (default is `FALSE`). If `TRUE`, the function will
  re-interpret the wildtype vs nonwildtype status for each observation
  based on the MIC values, against epidemiological cut-off (ECOFF)
  values. These will be reported in a new column `ecoff`, of class `sir`
  and coded as `NWT` (nonwildtype) or `WT` (wildtype).

## Value

Standardised AST data frame
