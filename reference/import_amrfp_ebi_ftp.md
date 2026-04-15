# Import EBI-processed AMRFinderPlus Genotypes from FTP

This function imports processed EBI-processed AMRFinderPlus genotyping
results. The expected input is genotype data retrieved from the [EBI AMR
Portal FTP
site](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) either
directly or via the function
[`download_ebi()`](https://amrgen.org/reference/download_ebi.md). Note
that files downloaded from the [EBI AMR Portal web
browser](https://www.ebi.ac.uk/amr/data/?view=predictions) are formatted
differently and can be imported using
[import_amrfp_ebi_web](https://amrgen.org/reference/import_amrfp_ebi_web.md).

## Usage

``` r
import_amrfp_ebi_ftp(input_table)
```

## Arguments

- input_table:

  R object or file path for the input EBI genotype table (R object, or
  file path to a TSV or CSV file).

## Value

A data frame with the processed genotype data, with harmonised gene
names, mapped drug agents, and drug classes which can be used for other
functions of the ARMgen package. The output retains the original columns
from the AMRFinderPlus table along with the newly mapped variables:

- `id`: The sample identifier (`character`).

- `marker`: The name of the genotype marker as it appears in the input
  (e.g. `gyrA_S83F`) (`character`).

- `gene`: The gene identifier (`character`).

- `mutation`: The mutation detected within the gene, converted to [HGVS
  nomenclature](https://hgvs-nomenclature.org/stable/) syntax (e.g.
  `Ser83Phe`) (`character`).

- `drug_class`: Name of the antibiotic group associated with the
  genotype marker, compatible with AMR pkg (`character`).

- `drug`: Name of the specific antibiotic associated with the genotype
  marker, compatible with AMR pkg (`ab`). Value `NA` is assigned when
  the markers are annotated with a class only and not a specific
  antibiotic. ... Other fields specific to the input file

## Details

These data are pre-processed by EBI to match NCBI class/subclass to
CARD's antibiotic resistance ontology (ARO), however for consistency
this function will re-process the data to generate `drug` and
`drug_class` fields consistent with the
[`import_amrfp()`](https://amrgen.org/reference/import_amrfp.md)
function (the EBI fields `antibiotic*` are also retained). Note several
AMRFinderPlus fields are excluded from EBI files, including hierarchy
node, method, percent identity and coverage; therefore unlike the
[`import_amrfp()`](https://amrgen.org/reference/import_amrfp.md)
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
ebi_geno_raw <- download_ebi(
  data = "genotype", species = "Escherichia coli",
  geno_subclass = "QUINOLONE"
)

# Format the file for import
ebi_geno <- import_amrfp_ebi_ftp(ebi_geno_raw)
} # }
```
