# Download NCBI antimicrobial susceptibility testing (AST) data

This function downloads antimicrobial susceptibility phenotype data from
the NCBI Pathogen Detection database via the BioSample API. Data are
retrieved in batches, parsed from XML, and returned as a tidy tibble
with metadata including BioSample ID, Bioproject ID, and organism name.

## Usage

``` r
download_ncbi_pheno(
  species,
  pheno_drug = NULL,
  max_records = 15000,
  batch_size = 200,
  sleep_time = 0.34,
  force_drug_name = FALSE,
  reformat = FALSE,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE
)
```

## Arguments

- species:

  Character. Organism name for the search query (e.g.,
  `"Salmonella enterica"`). Required.

- pheno_drug:

  Character or vector. Optional drug name/s to filter the returned data.
  Strings will be processed using the AMR package to standardize names
  before matching, so e.g. `"amikacin"` or `"Amikacin"` or `"ami"` will
  be parsed to "amikacin" before matching. This can be turned off by
  setting `force_drug_name=TRUE`. Full list of allowed drug names in
  NCBI: <https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/>.

- max_records:

  Integer. Maximum number of BioSample records to retrieve. Default is
  `15000`.

- batch_size:

  Integer. Number of records fetched per API request. Default is `200`
  which is recommended by NCBI.

- sleep_time:

  Numeric. Seconds to pause between batch requests to avoid overloading
  NCBI servers. Default is `0.34`.

- force_drug_name:

  Logical. If `TRUE`, turns off standardizing the drug name using the
  AMR package before filtering, so that matching is done exactly on the
  input string/s. Default is `FALSE`.

- reformat:

  Logical. If `TRUE`, reformats the output using
  [`import_ncbi_biosample()`](https://amrgen.org/reference/import_ncbi_biosample.md)
  for compatibility with AMR analysis workflows. Default is `FALSE`.
  When set to `TRUE`, the data can also be interpreted against
  breakpoints/ECOFF by setting the `interpret_*=TRUE`.

- interpret_eucast:

  Logical. Passed to
  [`import_ncbi_biosample()`](https://amrgen.org/reference/import_ncbi_biosample.md).
  If `TRUE`, interprets MIC values using EUCAST breakpoints. Default is
  `FALSE`. Only used if `reformat`=`TRUE`.

- interpret_clsi:

  Logical. Passed to
  [`import_ncbi_biosample()`](https://amrgen.org/reference/import_ncbi_biosample.md).
  If `TRUE`, interprets MIC values using CLSI breakpoints. Default is
  `FALSE`. Only used if `reformat`=`TRUE`.

- interpret_ecoff:

  Logical. Passed to
  [`import_ncbi_biosample()`](https://amrgen.org/reference/import_ncbi_biosample.md).
  If `TRUE`, interprets MIC values using ECOFF cutoffs. Default is
  `FALSE`. Only used if `reformat`=`TRUE`.

## Value

A tibble with one row per phenotype measure, with corresponding
BioSample metadata.

## Details

The function constructs an Entrez query of the form:
`"<organism> AND antibiogram[filter]"`. XML records are downloaded in
batches, parsed, and combined into a single table. The resulting tibble
contains phenotype test results and associated metadata including:

- `id`: BioSample identifier

- `BioProject`: BioProject accession ID

- `organism`: Organism name

- `Antibiotic`, `Phenotype`, `Measurement`, `Units`, `Method`, `System`,
  `Manufacturer`, `Panel`, `Standard`: phenotype data columns

The function can optionally filter by one or more drugs. It can also
optionally reformat data for compatibility with AMRgen functions via
[`import_ncbi_biosample()`](https://amrgen.org/reference/import_ncbi_biosample.md),
and interpret the raw data measures against breakpoints or ECOFF. See
[`import_ncbi_biosample()`](https://amrgen.org/reference/import_ncbi_biosample.md)
for details of output formats when these options are used.

## NCBI API usage

Users are encouraged to set an NCBI API key via
[`rentrez::set_entrez_key()`](https://docs.ropensci.org/rentrez/reference/set_entrez_key.html)
to increase request limits and comply with NCBI usage policies.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download phenotype data for Klebsiella quasipneumoniae
pheno <- download_ncbi_pheno("Klebsiella quasipneumoniae")

# Download Klebsiella quasipneumoniae data, filter to amikacin and ampicillin
pheno <- download_ncbi_pheno(
  "Klebsiella quasipneumoniae",
  pheno_drug = c("amikacin", "Amp")
)

# Download and reformat for AMRgen workflow with EUCAST interpretation
pheno <- download_ncbi_pheno(
  "Klebsiella quasipneumoniae",
  reformat = TRUE,
  interpret_eucast = TRUE
)
} # }
```
