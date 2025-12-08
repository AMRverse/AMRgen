# Get Binary Matrix of Genotype and Phenotype Data

This function generates a binary matrix representing the resistance (R
vs S/I) and nonwildtype (R/I vs S) status for a given antibiotic, and
presence or absence of genetic markers related to one or more specified
drug classes. It takes as input separate tables for genotype and
phenotype data, matches these according to a common identifier (either
specified by column names or assuming the first column contains the ID),
and filters the data according to the specified antibiotic and drug
class criteria before creating a binary matrix. Suitable input files can
be generated using `import_ncbi_ast` to import phenotype data from NCBI,
and `parse_amrfp` to import genotype data from AMRFinderPlus.

## Usage

``` r
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic,
  drug_class_list,
  keep_SIR = TRUE,
  keep_assay_values = FALSE,
  keep_assay_values_from = c("mic", "disk"),
  geno_sample_col = NULL,
  pheno_sample_col = NULL,
  sir_col = "pheno_clsi",
  ecoff_col = "ecoff",
  marker_col = "marker",
  most_resistant = TRUE
)
```

## Arguments

- geno_table:

  A data frame containing genotype data, including at least one column
  labeled `drug_class` for drug class information and one column for
  sample identifiers (specified via `geno_sample_col`, otherwise it is
  assumed the first column contains identifiers).

- pheno_table:

  A data frame containing phenotype data, which must include a column
  `drug_agent` (with the antibiotic information), a column with the
  resistance interpretation (S/I/R, specified via `sir_col`), and
  optionally a column with the ECOFF interpretation (WT/NWT, specified
  via `ecoff_col`).

- antibiotic:

  A character string specifying the antibiotic of interest to filter
  phenotype data. The value must match one of the entries in the
  `drug_agent` column of `pheno_table`.

- drug_class_list:

  A character vector of drug classes to filter genotype data for markers
  related to the specified antibiotic. Markers in `geno_table` will be
  filtered based on whether their `drug_class` matches any value in this
  list.

- keep_SIR:

  A logical indicating whether to retain the full S/I/R phenotype column
  in the output. Defaults to `TRUE`.

- keep_assay_values:

  A logical indicating whether to include columns with the raw phenotype
  assay data. Assumes there are columns labelled "mic" and "disk"; these
  will be added to the output table if present. Defaults to `FALSE`.

- keep_assay_values_from:

  A character vector specifying which assay values (e.g., `"mic"`,
  `"disk"`) to retain if `keep_assay_values` is `TRUE`. Defaults to
  `c("mic", "disk")`.

- geno_sample_col:

  A character string (optional) specifying the column name in
  `geno_table` containing sample identifiers. Defaults to `NULL`, in
  which case it is assumed the first column contains identifiers.

- pheno_sample_col:

  A character string (optional) specifying the column name in
  `pheno_table` containing sample identifiers. Defaults to `NULL`, in
  which case it is assumed the first column contains identifiers.

- sir_col:

  A character string specifying the column name in `pheno_table` that
  contains the resistance interpretation (SIR) data. The values should
  be interpretable as "R" (resistant), "I" (intermediate), or "S"
  (susceptible).

- ecoff_col:

  A character string specifying the column name in `pheno_table` that
  contains the ECOFF interpretation of phenotype. The values should be
  interpretable as "WT" (wildtype) or "NWT" (nonwildtype).

- marker_col:

  A character string specifying the column name in `geno_table`
  containing the marker identifiers. Defaults to `"marker"`.

- most_resistant:

  A logical indicating whether, when multiple phenotype entries are
  present for the same sample and drug, the most resistant should be
  kept (otherwise the least resistant is kept). Default is `TRUE`.

## Value

A data frame where each row represents a sample, and each column
represents a genetic marker related to the specified antibiotic's drug
class. The binary values in the matrix indicate the presence (`1`) or
absence (`0`) of each marker for each sample, along with resistance
status columns for the specified antibiotic: `R` for resistant (defined
from `sir_col`, 1=R, 0=I/S) and `NWT` for nonwildtype (defined by
`ecoff_col` if provided: 1=NWT, 0=WT; otherwise defined from `sir_col`:
1=I/R, 0=S).

## Details

This function performs several steps:

- Verifies that the `pheno_table` contains a `drug_agent` column and
  converts it to class `ab` if necessary.

- Filters the `pheno_table` to retain data related to the specified
  antibiotic.

- Checks that the `geno_table` contains markers associated with the
  specified drug class(es).

- Matches sample identifiers between `geno_table` and `pheno_table`.

- Extracts and transforms the phenotype data into a binary format
  indicating resistance and NWT status.

- Constructs a binary matrix for genotype data, with each column
  representing a genetic marker.

- Returns a single matrix with sample identifiers plus binary variables
  for each phenotype and genetic marker.

## See also

`compare_geno_pheno_id`, `as.ab`, `as.sir`, `ab_name`

## Examples

``` r
if (FALSE) { # \dontrun{
geno_table <- parse_amrfp("testdata/Ecoli_AMRfinderplus_n50.tsv", "Name")
pheno_table <- import_ncbi_ast("testdata/Ecoli_AST_NCBI_n50.tsv")
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype"
)
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype",
  keep_assay_values = TRUE
)
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype",
  keep_assay_values = TRUE,
  keep_assay_values_from = "mic"
)
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype",
  keep_assay_values = TRUE,
  keep_assay_values_from = "MIC (mg/L)"
)
} # }
```
