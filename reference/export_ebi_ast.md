# Export EBI Antibiogram

Format AMRgen long-format AST data to a table with the fields required
for submission to EBI, and optionally generate JSON submission files
(one per BioSample). See
<https://www.ebi.ac.uk/amr/amr_submission_guide/>).

## Usage

``` r
export_ebi_ast(
  data,
  pheno_col = "pheno_provided",
  breakpoint_version,
  submission_account,
  domain = "self.ExampleDomain",
  output_dir = NULL
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
  `method`, `platform`.

- pheno_col:

  Character string naming the column that contains SIR interpretations
  (class `sir`). Default `"pheno_provided"`.

- breakpoint_version:

  Character string specifying the breakpoint version used for
  interpretation (e.g. `"EUCAST 2024"`).

- submission_account:

  Character string specifying the EBI Webin submission account
  identifier (e.g. `"Webin-###"`). If not provided, JSON output files
  will not be generated and the function will return the formated table
  only, which can be further updated and converted to submission-ready
  JSON later using
  [`format_ebi_json()`](https://AMRverse.github.io/AMRgen/reference/format_ebi_json.md).

- domain:

  Character string specifying the domain used in the submission metadata
  (e.g. `"self.ExampleDomain"`). If not provided, JSON output files will
  not be generated and the function will return the formated table only,
  which can be further updated and converted to submission-ready JSON
  later using
  [`format_ebi_json()`](https://AMRverse.github.io/AMRgen/reference/format_ebi_json.md).

- output_dir:

  Character string specifying the directory where JSON files should be
  written. If not provided, JSON output files will not be generated and
  the function will return the formated table only, which can be further
  updated and converted to submission-ready JSON later using
  [`format_ebi_json()`](https://AMRverse.github.io/AMRgen/reference/format_ebi_json.md).

## Value

Formatted data frame. When `output_dir` is provided, the AST data is
also written to individual JSON submission files, one per BioSample, in
the specified directory.

## Details

Antibiotic names are in Title Case with `"/"` separating combination
agents (EBI convention, e.g. `"Amoxicillin/clavulanic acid"`).

Species names are derived from the `spp_pheno` column via
[`AMR::mo_name()`](https://amr-for-r.org/reference/mo_property.html).

## Examples

``` r
# Return formatted data frame without writing files
ebi_df <- export_ebi_ast(staph_ast_ebi)
if (FALSE) { # \dontrun{
# Write out data for each BioSample to an individual JSON file for submission
ebi_df <- export_ebi_ast(staph_ast_ebi,
  breakpoint_version = "EUCAST 2015",
  submission_account = "Webin-###",
  domain = "self.ExampleDomain",
  output_dir = "/path/to/output/"
)
} # }
```
