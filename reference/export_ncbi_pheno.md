# Export NCBI BioSample Antibiogram

Convert AMRgen long-format AST data to an [NCBI BioSample
Antibiogram](https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/)
submission file.

## Usage

``` r
export_ncbi_pheno(
  data,
  file = NULL,
  overwrite = FALSE,
  pheno_col = "pheno_provided",
  guideline = NULL,
  vendor = NULL,
  version = NULL
)
```

## Arguments

- data:

  A data frame in AMRgen long format (e.g. output of
  [`import_pheno()`](https://amrgen.org/reference/import_pheno.md) or
  [`format_pheno()`](https://amrgen.org/reference/format_pheno.md)).
  Expected columns: `id`, `drug`, and at least one phenotype column (see
  `pheno_col`). Optional columns: `mic`, `disk`, `method`, `guideline`,
  `platform`.

- file:

  File path for the output file (must end in `.txt` or `.tsv`). If
  `NULL` (default), no file is written and the formatted data frame is
  returned visibly.

- overwrite:

  Logical; overwrite an existing file? Default `FALSE`.

- pheno_col:

  Character string naming the column that contains SIR interpretations
  (class `sir`). Default `"pheno_provided"`.

- guideline:

  Optional single value to record in `testing_standard` field in the
  output (default `NULL`, in which case `testing_standard` will be
  populated from the `guideline` field in the input file).

- vendor:

  Optional single value to record in `vendor` field in the output
  (default `NULL`).

- version:

  Optional single value to record in
  `laboratory_typing_method_version_or_reagent` field in the output
  (default `NULL`).

## Value

When `file` is provided, the formatted data frame is returned invisibly
and a tab-delimited UTF-8 file is written to `file`. When `file = NULL`,
the formatted data frame is returned visibly and no file is written.

## Details

When both `mic` and `disk` columns are present, MIC values are preferred
(more precise). Disk values are only used for rows where MIC is `NA`.

MIC strings (e.g. `"<=0.5"`, `">=32"`, `"4"`) are split into a sign
(`<=`, `>=`, `<`, `>`, or `=`) and a numeric value.

Antibiotic names are converted to lowercase with combination separators
replaced by `"-"` (NCBI convention, e.g.
`"amoxicillin-clavulanic acid"`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Return formatted data frame without writing a file
ncbi_df <- export_ncbi_pheno(ecoli_pheno)

# Write out the ecoli_pheno data to file in NCBI format
export_ncbi_pheno(ecoli_pheno, "Ec_NCBI.tsv")

# Download data from EBI, then write it out to file in NCBI format
ebi_kq <- download_ebi(
  data = "phenotype",
  species = "Klebsiella quasipneumoniae",
  reformat = T
)
export_ncbi_pheno(ebi_kq, "Kq_NCBI.tsv")
} # }
```
