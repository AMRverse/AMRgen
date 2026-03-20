# Query antimicrobial genotype (MicroBIGG-E) data from NCBI Pathogen Detection BigQuery

This function queries the `ncbi-pathogen-detect.pdbrowser.microbigge`
BigQuery table to retrieve genotype data. **Note:** This function only
returns genotypes for BioSamples that also have AST data.

## Usage

``` r
query_ncbi_bq_geno(
  taxgroup,
  geno_subclass = NULL,
  geno_class = NULL,
  project_id = NULL
)
```

## Arguments

- taxgroup:

  String specifying the organism group to filter on (e.g., "Pseudomonas
  aeruginosa"). See <https://www.ncbi.nlm.nih.gov/pathogens/organisms/>
  for a list. Required.

- geno_subclass:

  (Optional) String or vector of strings specifying AMR subclasses to
  filter on (e.g., "CARBAPENEM").

- geno_class:

  (Optional) String or vector of strings specifying AMR classes to
  filter on (e.g., "BETA-LACTAM"). Ignored if `geno_subclass` is
  specified.

- project_id:

  (Optional) Google Cloud Project ID to use for billing. If NULL
  (default), looks for `GOOGLE_CLOUD_PROJECT` environment variable.

## Value

A tibble containing genotype data with columns renamed to match
[`import_amrfp()`](https://amrgen.org/reference/import_amrfp.md)
expectations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Query genotype data for Klebsiella pneumoniae samples that have AST data
geno_raw <- query_ncbi_bq_geno(
  taxgroup = "Klebsiella pneumoniae"
)

# Filter for carbapenem resistance genes and point mutations
geno_amrfp <- query_ncbi_bq_geno(
  taxgroup = "Klebsiella pneumoniae",
  geno_subclass = "CARBAPENEM"
)

geno <- import_amrfp(geno_amrfp, sample_col = "biosample_acc")
} # }
```
