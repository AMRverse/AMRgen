# Export NCBI BioSample Antibiogram

Convert AMRgen long-format AST data to an [NCBI BioSample
Antibiogram](https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/)
submission file.

## Usage

``` r
export_ncbi_biosample(
  data,
  file,
  overwrite = FALSE,
  pheno_col = "pheno_provided"
)
```

## Arguments

- data:

  A data frame in AMRgen long format (e.g. output of
  [`import_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ncbi_ast.md),
  [`import_ebi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast.md),
  or
  [`format_ast()`](https://AMRverse.github.io/AMRgen/reference/format_ast.md)).
  Expected columns: `id`, `drug_agent`, and at least one phenotype
  column (see `pheno_col`). Optional columns: `mic`, `disk`, `method`,
  `guideline`, `platform`.

- file:

  File path for the output file (must end in `.txt` or `.tsv`).

- overwrite:

  Logical; overwrite an existing file? Default `FALSE`.

- pheno_col:

  Character string naming the column that contains SIR interpretations
  (class `sir`). Default `"pheno_provided"`.

## Value

The formatted data frame is returned invisibly. A tab-delimited UTF-8
text file is written to `file`.

## Details

When both `mic` and `disk` columns are present, MIC values are preferred
(more precise). Disk values are only used for rows where MIC is `NA`.

MIC strings (e.g. `"<=0.5"`, `">=32"`, `"4"`) are split into a sign
(`<=`, `>=`, `<`, `>`, or `=`) and a numeric value.

Antibiotic names are converted to lowercase with combination separators
replaced by `"-"` (NCBI convention, e.g.
`"amoxicillin-clavulanic acid"`).
