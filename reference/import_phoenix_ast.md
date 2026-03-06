# Import and process antimicrobial phenotype data exported from BD Phoenix instruments

This function imports antimicrobial susceptibility testing (AST) data
from BD Phoenix instrument export files and converts them to the
standardised long-format used by AMRgen.

## Usage

``` r
import_phoenix_ast(
  input,
  sample_col = NULL,
  drug_col = NULL,
  mic_col = NULL,
  sir_col = NULL,
  species_col = NULL,
  species = NULL,
  ab = NULL,
  use_expertized = TRUE,
  source = NULL,
  instrument_guideline = NULL,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE
)
```

## Arguments

- input:

  Path to a Phoenix export file (XLS, XLSX, TXT, or TSV), or a data
  frame.

- sample_col:

  Name or 1-based column index of the sample ID column. If `NULL`
  (default), auto-detected from common names (`"Sample"`, `"ID"`,
  `"Isolate"`). For single-isolate files with no sample column, the
  filename is used as the sample ID.

- drug_col:

  Name or 1-based column index of the drug/antimicrobial column. If
  `NULL` (default), auto-detected from names such as `"Antimicrobial"`,
  `"Drug"`, etc.

- mic_col:

  Name or 1-based column index of the MIC column. If `NULL` (default),
  auto-detected from names such as `"MIC"`, `"MIC or Concentration"`,
  etc.

- sir_col:

  Name or 1-based column index of the S/I/R interpretation column. If
  `NULL` (default), auto-detected based on `use_expertized`: prefers
  `"Final (SIR)"`, `"Expert (SIR)"` when `TRUE`; prefers `"Interp"` when
  `FALSE`.

- species_col:

  Name or 1-based column index of the species/organism column. If `NULL`
  (default), auto-detected from names such as `"Organism"`, `"Species"`,
  etc. Ignored if `species` is supplied.

- species:

  Optional fixed species string applied to all rows (parsed via
  [`AMR::as.mo()`](https://amr-for-r.org/reference/as.mo.html)).
  Overrides `species_col`. Required for files without an organism
  column.

- ab:

  Optional antibiotic filter/override passed to
  [`interpret_ast()`](https://AMRverse.github.io/AMRgen/reference/interpret_ast.md).

- use_expertized:

  If `TRUE` (default), prefer the expert/final interpretation column
  over raw instrument interpretation when both are present and `sir_col`
  is not specified.

- source:

  Optional source label recorded in the `source` output column.

- instrument_guideline:

  Optional guideline string recorded in the `guideline` column (e.g.
  `"CLSI 2018"`, `"EUCAST 2024"`).

- interpret_eucast:

  Interpret MIC values against EUCAST human breakpoints (default
  `FALSE`).

- interpret_clsi:

  Interpret MIC values against CLSI human breakpoints (default `FALSE`).

- interpret_ecoff:

  Interpret MIC values against ECOFF values (default `FALSE`).

## Value

Standardised long-format AST data frame

## Details

Phoenix instruments can export data in different file formats. This
function handles both named-header exports (e.g. tab-delimited CLSI
reports with columns such as `"Antimicrobial"`,
`"MIC or Concentration"`, `"Final (SIR)"`) and positional/headerless
exports (XLS files with no column headers). Column detection is
automatic based on common Phoenix header name patterns, but any column
can be overridden by supplying its name or 1-based index via `drug_col`,
`mic_col`, `sir_col`, `sample_col`, and `species_col`.

For headerless files (detected when no standard column names are
recognised), the function assumes the standard Phoenix positional
layout: column 1 = sample ID, column 2 = species, column 3 = drug,
column 4 = MIC, column 5 = instrument interpretation, column 6 =
expert/final interpretation. Drug names are standardised using
[`AMR::as.ab()`](https://amr-for-r.org/reference/as.ab.html), which
handles abbreviated, full, and language-specific names (e.g.
`"Fosfomycin mit G6P"`). Species names are standardised using
[`AMR::as.mo()`](https://amr-for-r.org/reference/as.mo.html).

For tab-delimited files, trailing non-data sections commonly appended by
Phoenix (e.g. `"Resistance Markers"`, `"Expert Triggered Rules"`) are
automatically discarded.

## Examples

``` r
if (FALSE) { # \dontrun{
# Named-header export (e.g. CLSI report, auto-detected)
pheno <- import_phoenix_ast("phoenix_clsi_report.txt",
  species = "Escherichia coli", instrument_guideline = "CLSI 2018"
)

# Positional/headerless XLS export (auto-detected)
pheno <- import_phoenix_ast("phoenix_export.xls")

# Override column detection when names differ from defaults
pheno <- import_phoenix_ast("phoenix_export.xlsx",
  drug_col = "Antibiotic", mic_col = "Result",
  sir_col = "SIR_final", sample_col = "IsolateID",
  species_col = "Organism"
)
} # }
```
