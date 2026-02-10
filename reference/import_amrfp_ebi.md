# Import EBI-processed AMRFinderPlus Genotypes

This function imports EBI-processed AMRFinderPlus genotyping results.
The expected input is genotype data downloaded from the [EBI AMR Portal
web browser](https://www.ebi.ac.uk/amr/data/?view=predictions), or the
[EBI AMR Portal FTP
site](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) either
directly or via the function
[`download_ebi()`](https://AMRverse.github.io/AMRgen/reference/download_ebi.md).

## Usage

``` r
import_amrfp_ebi(input_table, web = FALSE)
```

## Arguments

- input_table:

  R object or file path for the input EBI genotype table (R object, or
  file path to a TSV or CSV file).

- web:

  Logical indicating whether the data is from the web portal (default
  `FALSE`). If `FALSE` input is assumed to be from FTP or
  [download_ebi](https://AMRverse.github.io/AMRgen/reference/download_ebi.md).

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
# Download quinolone-related genotype data for E. coli, from EBI
ebi_geno_raw <- download_ebi(data="genotype", species = "Escherichia coli", 
                        geno_subclass="QUINOLONE")

# Format the file for import
ebi_geno <- import_amrfp_ebi(ebi_geno_raw)

# Download data from EBI web portal and import the file
ebi_geno_from_web <- import_amrfp_ebi("amr_records.csv", web=TRUE)
} # }
```
