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
  the dataset (default `Name`).

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
  (`AFP_Subclass`) to standardised drug agents (`drug_agent`) and drug
  classes (`drug_class`). Defaults to `amrfp_drugs_table`, which is
  provided internally.

## Value

A tibble containing the processed AMR elements, with harmonised gene
names, mapped drug agents, and drug classes. The output retains the
original columns from the AMRFinderPlus table along with the newly
mapped variables.

## Details

The function performs the following steps:

- Reads the AMRFinderPlus output table.

- Filters the data to only include AMR elements.

- Converts gene symbols to a harmonised format.

- Splits multiple subclass annotations into separate rows.

- Maps AMRFinderPlus subclasses to standardised drug agent and drug
  class names using `amrfp_drugs`.

- Converts drug agent names to the `"ab"` class from the AMR package.
  This processing ensures compatibility with downstream AMR analysis
  workflows.

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
