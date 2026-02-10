# Import EBI-processed AMRFinderPlus Genotypes from Web

This function imports EBI-processed AMRFinderPlus genotyping results.
The expected input is genotype data downloaded from the [EBI AMR Portal
web browser](https://www.ebi.ac.uk/amr/data/?view=predictions). Note
that files downloaded from the [EBI AMR Portal FTP
site](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/), either
directly or via the function
[`download_ebi()`](https://AMRverse.github.io/AMRgen/reference/download_ebi.md),
are formatted differently and can be imported using
[import_amrfp_ebi_ftp](https://AMRverse.github.io/AMRgen/reference/import_amrfp_ebi_ftp.md).

## Usage

``` r
import_amrfp_ebi_web(input_table)
```

## Arguments

- input_table:

  R object or file path for the input EBI genotype table (R object, or
  file path to a TSV or CSV file).

## Value

A tibble containing the processed AMR elements, with harmonised gene
names, drug agents, and drug classes. The output retains the original
columns from the AMRFinderPlus table along with the newly mapped
variables.

## Details

These data are pre-processed by EBI to match NCBI class/subclass to
CARD's antibiotic resistance ontology (ARO), however for consistency
this function will re-process the data to generate `drug_agent` and
`drug_class` fields consistent with the
[`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
function (the EBI fields `antibiotic*` are also retained). Note several
AMRfinderplus fields are excluded from EBI files, including hierarchy
node, method, percent identity and coverage; therefore unlike the
[`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
function, this function cannot assign `variation type` or `node`.

The function performs the following steps:

- Reads the EBI-processed genotype table.

- Maps AMRFinderPlus subclasses to standardised drug agent and drug
  class names using `amrfp_drugs` (EBI-mappings are retained in
  `antibiotic*` fields.)

- Converts drug agent names to the `ab` class from the AMR package. This
  processing ensures compatibility with downstream AMR analysis
  workflows.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download data from EBI web portal and import the file
ebi_geno <- import_amrfp_ebi_web("amr_records.csv")
} # }
```
