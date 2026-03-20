# Query antimicrobial phenotype (AST) data from NCBI Pathogen Detection BigQuery

This function queries the `ncbi-pathogen-detect.pdbrowser.ast` BigQuery
table to retrieve antimicrobial susceptibility testing (AST) data from
NCBI Pathogen Detection.

## Usage

``` r
query_ncbi_bq_ast(
  taxgroup,
  antibiotic = NULL,
  force_antibiotic = FALSE,
  project_id = NULL
)
```

## Arguments

- taxgroup:

  String specifying the organism group to filter on (e.g., "Pseudomonas
  aeruginosa"). See <https://www.ncbi.nlm.nih.gov/pathogens/organisms/>
  for a list. Required.

- antibiotic:

  (Optional) String (or vector of strings) specifying the antibiotic
  name/s to filter on (default NULL). Uses the AMR package to try to fix
  typos, and format to lower-case.

- force_antibiotic:

  (Optional) Logical indicating whether to turn off parsing of
  antibiotic names and match exactly on the input strings (default
  FALSE).

- project_id:

  (Optional) Google Cloud Project ID to use for billing. If NULL
  (default), looks for `GOOGLE_CLOUD_PROJECT` environment variable.

## Value

A tibble containing AST data with columns renamed to match
[`import_ncbi_ast()`](https://amrgen.org/reference/import_ncbi_ast.md)
expectations.

## Details

Requires Google Cloud authentication. Run
[`bigrquery::bq_auth()`](https://bigrquery.r-dbi.org/reference/bq_auth.html)
before first use, or set up application default credentials via
`gcloud auth application-default login`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Query AST data for Klebsiella pneumoniae, filtered to meropenem
ast_raw <- query_ncbi_bq_ast(
  taxgroup = "Klebsiella pneumoniae",
  antibiotic = "meropenem"
)

# Import and reinterpret using CLSI breakpoints
ast <- import_ncbi_ast(ast_raw, interpret_clsi = TRUE)
} # }
```
