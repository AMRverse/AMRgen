# Export AST Data

Generic dispatcher that exports AMRgen long-format AST data to a
submission-ready file. Currently supports NCBI BioSample Antibiogram and
EBI Antibiogram formats.

## Usage

``` r
export_ast(
  data,
  file,
  format = "ncbi",
  overwrite = FALSE,
  pheno_col = "pheno_provided",
  ...
)
```

## Arguments

- data:

  A data frame in AMRgen long format.

- file:

  File path for the output file.

- format:

  Target format: `"ncbi"` (default) or `"ebi"`.

- overwrite:

  Logical; overwrite an existing file? Default `FALSE`.

- pheno_col:

  Character string naming the column that contains SIR interpretations.
  Default `"pheno_provided"`.

- ...:

  Additional arguments passed to the format-specific export function.

## Value

The formatted data frame is returned invisibly.
