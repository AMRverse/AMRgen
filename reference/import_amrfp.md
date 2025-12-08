# Import and Process AMRFinderPlus Results

This function imports and processes AMRFinderPlus results, extracting
antimicrobial resistance (AMR) elements and mapping them to standardised
antibiotic names and drug classes. The function also converts gene
symbols to a harmonised format and ensures compatibility with the AMR
package.

## Usage

``` r
import_amrfp(input_table, sample_col = "Name", amrfp_drugs = amrfp_drugs_table)
```

## Arguments

- input_table:

  A character string specifying the path to the AMRFinderPlus results
  table (TSV format).

- sample_col:

  A character string specifying the column that identifies samples in
  the dataset (default `Name`).

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
ecoli_geno_raw

# import first few rows of this data frame and parse it as AMRfp data
geno <- import_amrfp(ecoli_geno_raw %>% head(n = 10), "Name")
geno
} # }
```
