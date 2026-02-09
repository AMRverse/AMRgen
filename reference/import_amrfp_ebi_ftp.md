# Import EBI-processed AMRFinderPlus Genotypes

This function imports processed EBI-processed AMRFinderPlus genotyping
results. The expected input is genotype data retrieved from the [EBI AMR
Portal FTP site](ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/)
either directly or via the function
[`download_ebi()`](https://AMRverse.github.io/AMRgen/reference/download_ebi.md).
Note that files downloaded from the [EBI AMR Portal web
browser](https://www.ebi.ac.uk/amr/data/?view=predictions) are formatted
differently.

## Usage

``` r
import_amrfp_ebi_ftp(
  input_table,
  sample_col = "BioSample_ID",
  element_symbol_col = "amr_element_symbol",
  element_type_col = "element_type",
  element_subtype_col = "element_subtype",
  subclass_col = "subclass",
  class_col = "class",
  amrfp_drugs = amrfp_drugs_table
)
```

## Arguments

- input_table:

  A character string specifying the path to the EBI genotype table (R
  object, or file path to a TSV file).

- sample_col:

  A character string specifying the column that identifies samples in
  the dataset (default "BioSample_ID").

- element_symbol_col:

  A character string specifying the column that contains the content of
  the AMRfinderplus 'Element symbol' field (default
  "amr_element_symbol").

- element_type_col:

  A character string specifying the column that contains the content of
  the AMRfinderplus 'Element type' field (default "element_type").

- element_subtype_col:

  A character string specifying the column that contains the content of
  the AMRfinderplus 'Element subtype' field (default "element_subtype").

- subclass_col:

  A character string specifying the column that contains the content of
  the AMRfinderplus 'Subclass' field (default "subclass").

- class_col:

  A character string specifying the column that contains the content of
  the AMRfinderplus 'Class' field (default "class").

- amrfp_drugs:

  Data from mapping drug names to classes (default amrfp_drugs_table).

## Value

A tibble containing the processed AMR elements, with harmonised gene
names, drug agents (mapped by ), and drug classes. The output retains
the original columns from the AMRFinderPlus table along with the newly
mapped variables.

## Details

These data are pre-processed by EBI to match NCBI class/subclass to
CARD's antibiotic resistance ontology (ARO), however for consistency
this function will re-process the data to generate 'drug_agent' and
'drug_class' fields consistent with the
[`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
function (the EBI fields 'antibiotic\*' are also retained). Note several
AMRfinderplus fields are excluded from EBI files, including hierarchy
node, method, percent identity and coverage; therefore unlike the
[`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
function, this function cannot assign 'variation type' or 'node'.

The function performs the following steps:

- Reads the EBI-processed genotype table.

- Maps AMRFinderPlus subclasses to standardised drug agent and drug
  class names using `amrfp_drugs` (EBI-mappings are retained in
  ('antibiotic\*' fields.)

- Converts drug agent names to the `"ab"` class from the AMR package.
  This processing ensures compatibility with downstream AMR analysis
  workflows.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download quinolone-related genotype data for E. coli, from EBI
ebi_geno_raw <- download_ebi(data="genotype", species = "Escherichia coli", 
                        geno_subclass="QUINOLONE")

# Format the file for import
ebi_geno <- import_amrfp_ebi_ftp(ebi_geno_raw)
} # }
```
