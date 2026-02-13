# Get Binary Matrix of Genotype and Phenotype Data

This function generates a binary matrix representing the resistance (R
vs S/I) and nonwildtype (NWT vs WT, or R/I vs S) status for a given
antibiotic, and presence or absence of genetic markers related to one or
more specified drug classes. It takes as input separate tables for
genotype and phenotype data, matches these according to a common
identifier (either specified by column names or assuming the first
column contains the ID), and filters the data according to the specified
antibiotic and drug class criteria before creating a binary matrix.
Suitable input files can be generated using
[`import_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ast.md)
to import phenotype data, and
[`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
to import genotype data from AMRFinderPlus.

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

  A data frame containing genotype data, in long form with one row per
  sample and genetic marker. Expected format is that output by
  [`import_amrfp()`](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md)
  and must include a column labeled `drug_class` (indicating the
  antibiotic class associated with each marker), in addition to a column
  indicating the marker (column name specified via `marker_col`) and a
  column for sample identifiers (specified via `geno_sample_col`,
  otherwise it is assumed the first column contains identifiers).

- pheno_table:

  A data frame containing phenotype data, in long form with one row per
  sample, drug and assay result. Expected format is that output by
  [`import_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ast.md)
  and must include a column `drug_agent` (indicating the drug agent,
  interpretable as AMR pkg class `ab`), in addition to a column for
  sample identifiers (specified via `pheno_sample_col`, otherwise it is
  assumed the first column contains identifiers), a column with the
  resistance interpretation (S/I/R, specified via `sir_col`), and
  optionally a column with the ECOFF interpretation (WT/NWT or S/R,
  specified via `ecoff_col`).

- antibiotic:

  A character string specifying the antibiotic of interest to filter
  phenotype data. The value must match one of the entries in the
  `drug_agent` column of `pheno_table` or be coercible to a match using
  [AMR::as.ab](https://amr-for-r.org/reference/as.ab.html).

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
  assay data. Assumes there are columns labelled `"mic"` and `"disk"`;
  these will be added to the output table if present. Defaults to
  `FALSE`.

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
  be interpretable as `"R"`, `"I"`, `"S"`. Default `"pheno_clsi`.

- ecoff_col:

  A character string specifying the column name in `pheno_table` that
  contains the ECOFF interpretation of phenotype. The values should be
  interpretable as `"WT"` (wildtype) and `"NWT"` (nonwildtype), or `"S"`
  / `"I"` / `"R"`. Default `"ecoff`.

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
# Import example E. coli AMRFinderPlus data from AllTheBacteria
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
geno_pheno_cip <- get_binary_matrix(
  ecoli_geno,
  ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "pheno_clsi"
)
geno_pheno_cip <- get_binary_matrix(
  ecoli_geno,
  ecoli_ast,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype",
  keep_assay_values = TRUE
)
} # }
```
