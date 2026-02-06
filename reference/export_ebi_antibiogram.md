# Export EBI Antibiogram

Convert AMRgen long-format AST data to an EBI antibiogram submission
file (see [EBI
COMPARE-AMR](https://github.com/EBI-COMMUNITY/compare-amr)).

## Usage

``` r
export_ebi_antibiogram(
  data,
  file,
  overwrite = FALSE,
  pheno_col = "pheno_provided",
  sep = "\t"
)
```

## Arguments

- data:

  A data frame in AMRgen long format (e.g. output of
  [`import_ebi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast.md),
  [`import_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ncbi_ast.md),
  or
  [`format_ast()`](https://AMRverse.github.io/AMRgen/reference/format_ast.md)).
  Expected columns: `id`, `drug_agent`, `spp_pheno`, and at least one
  phenotype column (see `pheno_col`). Optional columns: `mic`, `disk`,
  `method`, `guideline`, `platform`.

- file:

  File path for the output file.

- overwrite:

  Logical; overwrite an existing file? Default `FALSE`.

- pheno_col:

  Character string naming the column that contains SIR interpretations
  (class `sir`). Default `"pheno_provided"`.

- sep:

  Field separator for the output file. Default `"\t"` (tab-delimited).
  Use `","` for CSV.

## Value

The formatted data frame is returned invisibly. A file is written to
`file`.

## Details

Antibiotic names are in Title Case with `"/"` separating combination
agents (EBI convention, e.g. `"Amoxicillin/clavulanic acid"`).

Species names are derived from the `spp_pheno` column via
[`AMR::mo_name()`](https://amr-for-r.org/reference/mo_property.html).
