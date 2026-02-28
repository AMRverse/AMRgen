# Export AST Data

Generic dispatcher that exports AMRgen long-format AST data to a
submission-ready file. Currently supports NCBI BioSample Antibiogram and
EBI Antibiogram formats.

## Usage

``` r
export_ast(
  data,
  file = NULL,
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

  File path for the output file. If `NULL` (default), no file is written
  and the formatted data frame is returned visibly.

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

When `file` is provided, the formatted data frame is returned invisibly
and a file is written. When `file = NULL`, the formatted data frame is
returned visibly and no file is written.

## Examples

``` r
if (FALSE) { # \dontrun{
# Return NCBI formatted data frame without writing a file
ncbi_df <- export_ast(ecoli_ast)

# Write out the ecoli_ast data to file in EBI format
export_ast(ecoli_ast, "Ec_EBI.tsv", format = "ebi")

# Download data from EBI, then write it out to file in NCBI format
ebi_kq <- download_ebi(
  data = "phenotype",
  species = "Klebsiella quasipneumoniae",
  reformat = T
)
export_ast(ebi_kq, "Kq_NCBI.tsv", format = "ncbi")
} # }
```
