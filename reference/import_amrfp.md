# Import and Process AMRFinderPlus Results

This function imports and processes AMRFinderPlus results, extracting
antimicrobial resistance (AMR) elements and mapping them to standardised
antibiotic names and drug classes. The function also converts gene
symbols to a harmonised format and ensures compatibility with the AMR
package.

## Usage

``` r
import_amrfp(
  input_table,
  sample_col = "Name",
  element_symbol_col = NULL,
  element_type_col = NULL,
  element_subtype_col = "Element subtype",
  method_col = "Method",
  node_col = "Hierarchy node",
  subclass_col = "Subclass",
  class_col = "Class",
  amrfp_drugs = amrfp_drugs_table
)
```

## Arguments

- input_table:

  A character string specifying a dataframe or path to the AMRFinderPlus
  results table (TSV format).

- sample_col:

  A character string specifying the column that identifies samples in
  the dataset (default "`Name`").

- element_symbol_col:

  Optional character string specifying the column containing gene or
  element symbols if non-standard column names are used.

- element_type_col:

  Optional character string specifying the column indicating element
  type (e.g. AMR).

- element_subtype_col:

  Character string specifying the column used to detect mutation
  subtypes.

- method_col:

  Character string specifying the AMRFinderPlus method column.

- node_col:

  Character string specifying the hierarchy node column.

- subclass_col:

  Character string specifying the AMRFinderPlus subclass column.

- class_col:

  Character string specifying the AMRFinderPlus class column.

- amrfp_drugs:

  A tibble containing a reference table mapping AMRFinderPlus subclasses
  (`AMRFP_Subclass`) to standardised drug classes (`drug_class`).
  Defaults to `amrfp_drugs_table`, which is provided internally.

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

- `node`: The node in the NCBI reference gene hierarchy corresponding to
  the gene (`character`).

- `drug_class`: Name of the antibiotic group associated with the
  genotype marker, compatible with AMR pkg (`character`).

- `drug`: Name of the specific antibiotic associated with the genotype
  marker, compatible with AMR pkg (`ab`). Value `NA` is assigned when
  the markers are annotated with a class only and not a specific
  antibiotic.

- `variation type`: Type of variation, e.g. `Gene presence detected`,
  `Protein variant detected`, `Nucleotide variant detected`,
  `Inactivating mutation detected`, `Promoter variant detected`. ...
  Other fields specific to the input file

## Details

The function performs the following steps:

- Reads the AMRFinderPlus output table.

- Filters the data to only include AMR elements.

- Converts gene symbols to a harmonised format.

- Splits multiple subclass annotations into separate rows.

- Maps AMRFinderPlus subclasses to standardised drug class names
  recognised by the AMR pkg. This processing ensures compatibility with
  downstream AMRgen analysis workflows.

## Examples

``` r
if (FALSE) { # \dontrun{
# small example E. coli AMRFinderPlus data
data(ecoli_geno_raw)
ecoli_geno_raw

# import first few rows of this data frame and parse it as AMRfp data
geno <- import_amrfp(ecoli_geno_raw %>% head(n = 10), "Name")
geno
} # }
```
