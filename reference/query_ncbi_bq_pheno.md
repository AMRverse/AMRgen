# Query antimicrobial phenotype data from NCBI Pathogen Detection BigQuery

This function queries the `ncbi-pathogen-detect.pdbrowser.ast` BigQuery
table to retrieve antimicrobial phenotype data from NCBI Pathogen
Detection.

## Usage

``` r
query_ncbi_bq_pheno(
  taxgroup,
  pheno_drug = NULL,
  force_drug_name = FALSE,
  project_id = NULL
)
```

## Arguments

- taxgroup:

  String specifying the organism group to filter on (e.g., "Pseudomonas
  aeruginosa"). See <https://www.ncbi.nlm.nih.gov/pathogens/organisms/>
  for a list. Required.

- pheno_drug:

  (Optional) String (or vector of strings) specifying the drug name/s to
  filter on (default NULL). Uses the AMR package to try to fix typos,
  and format to lower-case.

- force_drug_name:

  (Optional) Logical indicating whether to turn off parsing of drug
  names and match exactly on the input strings (default `FALSE`).

- project_id:

  (Optional) Google Cloud Project ID to use for billing. If NULL
  (default), looks for `GOOGLE_CLOUD_PROJECT` environment variable.

## Value

A tibble containing phenotype data with columns renamed to match
[`import_ncbi_pheno()`](https://amrgen.org/reference/import_ncbi_pheno.md)
expectations.

## Details

Requires Google Cloud authentication. Run
[`bigrquery::bq_auth()`](https://bigrquery.r-dbi.org/reference/bq_auth.html)
before first use, or set up application default credentials via
`gcloud auth application-default login`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Query phenotype data for Klebsiella pneumoniae, filtered to meropenem
pheno_raw <- query_ncbi_bq_pheno(
  taxgroup = "Klebsiella pneumoniae",
  pheno_drug = "meropenem"
)

# Import and reinterpret using CLSI breakpoints
pheno <- import_ncbi_pheno(pheno_raw, interpret_clsi = TRUE)
} # }
```
