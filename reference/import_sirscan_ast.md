# Import and process antimicrobial susceptibility data from SIRscan (Bio-Rad)

SIRscan (Bio-Rad) exports AST results as three separate
semicolon-delimited CSV files: MIC values (CMI), disk diffusion
diameters (DIAM), and SIR interpretations (INTERPR). Each file shares
the same row structure (one row per isolate) and antibiotic columns.
This function reads one or more of these files, pivots to long format,
and returns a standardised AMRgen data frame.

## Usage

``` r
import_sirscan_ast(
  mic_file = NULL,
  disk_file = NULL,
  interpr_file = NULL,
  source = NULL,
  species = NULL,
  ab = NULL,
  instrument_guideline = NULL,
  sirscan_codes = sirscan_codes,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE
)
```

## Arguments

- mic_file:

  Path to the SIRscan CMI (MIC) export file. Optional.

- disk_file:

  Path to the SIRscan DIAM (disk diffusion) export file. Optional.

- interpr_file:

  Path to the SIRscan INTERPR (SIR interpretation) export file.
  Optional. At least one of `mic_file`, `disk_file`, or `interpr_file`
  must be supplied.

- source:

  Optional source label for the `source` output column.

- species:

  Optional fixed species string applied to all rows (parsed via
  [`AMR::as.mo()`](https://amr-for-r.org/reference/as.mo.html)). If
  `NULL`, the `Germe` column from the file is used.

- ab:

  Optional antibiotic filter/override passed to
  [`interpret_ast()`](https://amrgen.org/reference/interpret_ast.md).

- instrument_guideline:

  Optional guideline string for the `guideline` column (e.g.
  `"EUCAST 2023"`).

- sirscan_codes:

  A data frame mapping SIRscan antibiotic codes (`sirscan_code`) to
  antibiotic names (`ab_name`) recognised by
  [`AMR::as.ab()`](https://amr-for-r.org/reference/as.ab.html). Defaults
  to the built-in `sirscan_codes` table. Supply your own data frame to
  override mappings — for example if your SIRscan configuration uses
  `S`, `K`, or `C` for different agents than the defaults. Set to `NULL`
  to skip the lookup and pass all codes directly to
  [`AMR::as.ab()`](https://amr-for-r.org/reference/as.ab.html).

- interpret_eucast:

  Interpret MIC/disk values against EUCAST human breakpoints (default
  `FALSE`).

- interpret_clsi:

  Interpret MIC/disk values against CLSI human breakpoints (default
  `FALSE`).

- interpret_ecoff:

  Interpret MIC/disk values against ECOFF values (default `FALSE`).

## Value

A standardised long-format data frame with columns: `id`, `drug_agent`,
`mic`, `disk`, `pheno_provided`, `pheno_eucast`, `pheno_clsi`, `ecoff`,
`method`, `platform`, `guideline`, `source`, `spp_pheno`,
`collection_date`, `specimen_type`.

## Details

Files have a 4-row metadata header (institution, lab, export date,
period) followed by a semicolon-delimited data block. Latin-1
(Windows-1252) encoding is assumed, as is standard for French-locale
SIRscan exports.

## Examples

``` r
if (FALSE) { # \dontrun{
# Import all three SIRscan export files
pheno <- import_sirscan_ast(
  mic_file     = "CMI_SALMO.csv",
  disk_file    = "DIAM_SALMO.csv",
  interpr_file = "INTERPR_SALMO.csv"
)

# MIC + interpretations only, with EUCAST re-interpretation
pheno <- import_sirscan_ast(
  mic_file = "CMI_SALMO.csv",
  interpr_file = "INTERPR_SALMO.csv",
  species = "Salmonella enterica",
  interpret_eucast = TRUE
)
} # }
```
