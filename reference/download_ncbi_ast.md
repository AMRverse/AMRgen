# Download NCBI antimicrobial susceptibility testing (AST) data

This function downloads antimicrobial susceptibility testing (AST) data
from the NCBI Pathogen Detection database via the BioSample API. Data
are retrieved in batches, parsed from XML, and returned as a tidy tibble
with metadata including BioSample ID, Bioproject ID, and organism name.

## Usage

``` r
download_ncbi_ast(
  species,
  antibiotic = NULL,
  max_records = 15000,
  batch_size = 200,
  sleep_time = 0.34,
  force_antibiotic = FALSE,
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

- antibiotic:

  Character or vector. Optional antibiotic name/s to filter the returned
  data. Strings will be processed using the AMR package to standardize
  names before matching, so e.g. `"amikacin"` or `"Amikacin"` or `"ami"`
  will be parsed to "amikacin" before matching. This can be turned off
  by setting `force_antibiotic=TRUE`. Full list of allowed antibiotic
  names in NCBI:
  <https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/>.

- max_records:

  Integer. Maximum number of BioSample records to retrieve. Default is
  `15000`.

- batch_size:

  Integer. Number of records fetched per API request. Default is `200`
  which is recommended by NCBI.

- sleep_time:

  Numeric. Seconds to pause between batch requests to avoid overloading
  NCBI servers. Default is `0.34`.

- force_antibiotic:

  Logical. If `TRUE`, turns off standardizing the antibiotic name using
  the AMR package before filtering, so that matching is done exactly on
  the input string/s. Default is `FALSE`.

- reformat:

  Logical. If `TRUE`, reformats the output using
  [import_ncbi_biosample](https://AMRverse.github.io/AMRgen/reference/import_ncbi_biosample.md)
  for compatibility with AMR analysis workflows. Default is `FALSE`.
  When set to `TRUE`, the data can also be interpreted against
  breakpoints/ECOFF by setting the `interpret_*=TRUE`.

- interpret_eucast:

  Logical. Passed to
  [interpret_ast](https://AMRverse.github.io/AMRgen/reference/interpret_ast.md)
  via
  [import_ncbi_biosample](https://AMRverse.github.io/AMRgen/reference/import_ncbi_biosample.md).
  If `TRUE`, interprets MIC values using EUCAST breakpoints. Default is
  `FALSE`. Only used if `reformat`=`TRUE`.

- interpret_clsi:

  Logical. Passed to
  [interpret_ast](https://AMRverse.github.io/AMRgen/reference/interpret_ast.md)
  via
  [import_ncbi_biosample](https://AMRverse.github.io/AMRgen/reference/import_ncbi_biosample.md).
  If `TRUE`, interprets MIC values using CLSI breakpoints. Default is
  `FALSE`. Only used if `reformat`=`TRUE`.

- interpret_ecoff:

  Logical. Passed to
  [interpret_ast](https://AMRverse.github.io/AMRgen/reference/interpret_ast.md)
  via
  [import_ncbi_biosample](https://AMRverse.github.io/AMRgen/reference/import_ncbi_biosample.md).
  If `TRUE`, interprets MIC values using ECOFF cutoffs. Default is
  `FALSE`. Only used if `reformat`=`TRUE`.

## Value

A tibble with one row per AST measure, with corresponding BioSample
metadata.

## Details

The function constructs an Entrez query of the form:
`"<organism> AND antibiogram[filter]"`. XML records are downloaded in
batches, parsed, and combined into a single table. The resulting tibble
contains AST test results and associated metadata including:

- `id`: BioSample identifier

- `BioProject`: BioProject accession ID

- `organism`: Organism name

- `Antibiotic`, `Phenotype`, `Measurement`, `Units`, `Method`, `System`,
  `Manufacturer`, `Panel`, `Standard`: AST data columns

The function can optionally filter by a one or more antibiotics. It can
also optionally reformat data for compatibility with AMRgen functions
via
[import_ncbi_ast](https://AMRverse.github.io/AMRgen/reference/import_ncbi_ast.md),
and interpret the raw data measures against breakpoints or ECOFF. See
[import_ncbi_ast](https://AMRverse.github.io/AMRgen/reference/import_ncbi_ast.md)
for details of output formats when these options are used.

## NCBI API usage

Users are encouraged to set an NCBI API key via
[`rentrez::set_entrez_key()`](https://docs.ropensci.org/rentrez/reference/set_entrez_key.html)
to increase request limits and comply with NCBI usage policies.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download AST data for Klebsiella quasipneumoniae
ast <- download_ncbi_ast("Klebsiella quasipneumoniae")

# Download Klebsiella quasipneumoniae data, filter to amikacin and ampicillin
ast <- download_ncbi_ast(
  "Klebsiella quasipneumoniae",
  antibiotic = c("amikacin", "Amp")
)

# Download and reformat for AMRgen workflow with EUCAST interpretation
ast <- download_ncbi_ast(
  "Klebsiella quasipneumoniae",
  reformat = TRUE,
  interpret_eucast = TRUE
)
} # }
```
