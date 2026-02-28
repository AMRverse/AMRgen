# Export EBI Antibiogram

Convert AMRgen long-format AST data to an EBI antibiogram submission
file (see [EBI
COMPARE-AMR](https://github.com/EBI-COMMUNITY/compare-amr)).

## Usage

``` r
export_ebi_antibiogram(
  data,
  file = NULL,
  overwrite = FALSE,
  pheno_col = "pheno_provided",
  sep = "\t"
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

- overwrite:

  Logical; overwrite an existing file? Default `FALSE`.

- pheno_col:

  Character string naming the column that contains SIR interpretations
  (class `sir`). Default `"pheno_provided"`.

- sep:

  Field separator for the output file. Default `"\t"` (tab-delimited).
  Use `","` for CSV.

## Value

When `file` is provided, the formatted data frame is returned invisibly
and a file is written to `file`. When `file = NULL`, the formatted data
frame is returned visibly and no file is written.

## Details

Antibiotic names are in Title Case with `"/"` separating combination
agents (EBI convention, e.g. `"Amoxicillin/clavulanic acid"`).

Species names are derived from the `spp_pheno` column via
[`AMR::mo_name()`](https://amr-for-r.org/reference/mo_property.html).

## Examples

``` r
if (FALSE) { # \dontrun{
# Return formatted data frame without writing a file
ebi_df <- export_ebi_antibiogram(ecoli_ast)

# Write out the ecoli_ast data to file in EBI format
export_ebi_antibiogram(ecoli_ast, "Ec_EBI.tsv")
} # }
```
