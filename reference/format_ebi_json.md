# Generate EBI antibiogram submission in JSON

Converts the tabular output of
[export_ebi_antibiogram](https://AMRverse.github.io/AMRgen/reference/export_ebi_antibiogram.md)
into JSON files formatted for submission to EBI as BioSample data
(https://www.ebi.ac.uk/amr/amr_submission_guide/). Each row of the input
dataset is converted into JSON records and printed to file.

## Usage

``` r
format_ebi_json(
  ebi_antibiogram_table,
  breakpoint_version,
  submission_account,
  domain = "self.ExampleDomain",
  output_dir
)
```

## Arguments

- ebi_antibiogram_table:

  A data frame in the format output by
  [export_ebi_antibiogram](https://AMRverse.github.io/AMRgen/reference/export_ebi_antibiogram.md).

- breakpoint_version:

  Character string specifying the breakpoint version used for
  interpretation (e.g. `"EUCAST 2024"`).

- submission_account:

  Character string specifying the Webin submission account identifier
  (e.g. `"Webin-###"`).

- domain:

  Character string specifying the domain used in the submission metadata
  (default `"self.ExampleDomain"`).

- output_dir:

  Character string specifying the directory where JSON files should be
  written.

## Value

Invisibly returns `NULL`. The function prints JSON-formatted AMR
submission records to file.

## Details

The function iterates over each biosample in `ebi_antibiogram_table` and
constructs a nested JSON object describing the antimicrobial
susceptibility testing result. Each record contains antibiotic metadata,
AST standards, measurement values, and resistance phenotype information.

JSON formatting is performed using
[`jsonlite::toJSON()`](https://jeroen.r-universe.dev/jsonlite/reference/fromJSON.html)
with `pretty = TRUE` and `auto_unbox = TRUE` .

## Examples

``` r
if (FALSE) { # \dontrun{
format_ebi_json(
  ast_dataset,
  breakpoint_version = "EUCAST 2015",
  submission_account = "Webin-###",
  domain = "self.ExampleDomain",
  output_dir = "/path/to/output/"
)
} # }
```
