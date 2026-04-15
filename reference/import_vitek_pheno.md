# Import and process antimicrobial phenotype data exported from Vitek instruments

This function imports antimicrobial susceptibility testing (AST) data
from Vitek instrument output files (wide CSV format) and converts it to
the standardised long-format used by AMRgen.

## Usage

``` r
import_vitek_pheno(
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

  String indicating the name of the column containing sample
  identifiers. Default: `"Lab ID"`

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

- use_expertized:

  Use expertized SIR (`TRUE`, default) or instrument SIR (`FALSE`)

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

- include_dates:

  Logical indicating whether to include `collection_date` and
  `testing_date` in output. Default: `TRUE`

## Value

Standardised AST data frame
