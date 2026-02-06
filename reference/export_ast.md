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

  A data frame in AMRgen long format (e.g. output of
  [`import_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ast.md)
  or
  [`format_ast()`](https://AMRverse.github.io/AMRgen/reference/format_ast.md)).
  Expected columns: `id`, `drug_agent`, `spp_pheno`, and at least one
  phenotype column (see `pheno_col`). Optional columns: `mic`, `disk`,
  `method`, `guideline`, `platform`.

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
