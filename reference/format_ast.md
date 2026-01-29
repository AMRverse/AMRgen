# Import and Process AST Data from a generic format

This function attempts to import antibiotic susceptibility testing (AST)
data in long-form antibiogram format (one row per sample and test),
suitable for downstream use with AMRgen analysis functions. It assumes
that the input file is a tab-delimited text file (e.g., TSV) or CSV
(which may be compressed) and parses relevant columns (antibiotic names,
species names, MIC or disk data, S/I/R calls) into suitable classes
using the AMR package. It optionally can use the AMR package to
interpret susceptibility phenotype (SIR) based on EUCAST or CLSI
guidelines (human breakpoints and/or ECOFF). If expected columns are not
found warnings will be given, and interpretation may not be possible.

## Usage

``` r
format_ast(
  input,
  sample_col = "id",
  species = NULL,
  species_col = "spp_pheno",
  ab = NULL,
  ab_col = "drug_agent",
  mic_col = "mic",
  disk_col = "disk",
  pheno_cols = c("ecoff", "pheno_eucast", "pheno_clsi", "pheno_provided"),
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE,
  rename_cols = TRUE
)
```

## Arguments

- input:

  A string representing a dataframe, or a path to an input file,
  containing the AST data in long-form antibiogram format (one row per
  sample and test). This might be a file containing the content of an
  EBI/NCBI AST dataset previously processed using
  [`import_ebi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast.md),
  [`import_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ncbi_ast.md),
  or
  [`import_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ast.md)
  functions, or files with a similar format/structure but with different
  column names.

- sample_col:

  (optional, default 'id') Name of the input data column that provides
  the sample name. If the 'rename' parameter is set to TRUE, this column
  will be renamed as 'id'.

- species:

  (optional) Name of the single species to which all samples belong. Use
  this if you want to interpret assay measurements but the input file
  does not contain a column indicating the species for each sample
  (called 'species' or a name specified by the `species_col` parameter).

- species_col:

  (optional, default 'species') Name of the input data column that
  provides a species name. If provided, this column will be converted to
  micro-organism class 'mo' via `as.mo()`. If the 'rename' parameter is
  set to TRUE, this column will also be renamed as 'spp_pheno'. If
  interpretation is switched on, this column will be used to identify
  the appropriate breakpoints for interpretation of each row in the data
  table.

- ab:

  (optional) Name of a single antibiotic to use for phenotype
  interpretation. Use this if you want to interpret assay measurements
  but the input file does not contain a column indicating the drug for
  each sample (called 'drug_agent' or a name specified by the `ab_col`
  parameter).

- ab_col:

  (optional, default 'drug_agent') Name of the input data column that
  provides a drug name. If provided, this column will be converted to
  antibiotic class 'ab' via `as.ab()`. If the 'rename' parameter is set
  to TRUE, this column will also be renamed as 'drug_agent'. If
  interpretation is switched on, this column will be used to identify
  the appropriate breakpoints for interpretation of each row in the data
  table.

- mic_col:

  (optional, default 'mic') Name of the input data column that provides
  MIC measurements. If provided, this column will be converted to MIC
  class 'mic' via `as.mic()`. If the 'rename' parameter is set to TRUE,
  this column will also be renamed as 'mic'. If interpretation is
  switched on, the MIC values will be interpreted against clinical
  breakpoints.

- disk_col:

  (optional, default 'disk') Name of the input data column that provides
  disk diffusion zone measurements. If provided, this column will be
  converted to disk diffusion class 'disk' via `as.disk()`. If the
  'rename' parameter is set to TRUE, this column will also be renamed as
  'disk'. If interpretation is switched on, the zone values will be
  interpreted against clinical breakpoints.

- pheno_cols:

  (optional, default
  `c("ecoff", "pheno_eucast", "pheno_clsi", "pheno_provided")`) Name of
  the input data column/s that provides disk diffusion zone measurements
  (as a character vector, or single string for a single column). If
  provided, these columns will be converted to SIR class 'sir' via
  `as.sir()`.

- interpret_eucast:

  A logical value (default is FALSE). If `TRUE`, the function will
  interpret the susceptibility phenotype (SIR) for each row based on the
  MIC or disk diffusion values, against EUCAST human breakpoints. These
  will be reported in a new column `pheno_eucast`, of class 'sir'.

- interpret_clsi:

  A logical value (default is FALSE). If `TRUE`, the function will
  interpret the susceptibility phenotype (SIR) for each row based on the
  MIC or disk diffusion values, against CLSI human breakpoints. These
  will be reported in a new column `pheno_clsi`, of class 'sir'.

- interpret_ecoff:

  A logical value (default is FALSE). If `TRUE`, the function will
  interpret the wildtype vs nonwildtype status for each row based on the
  MIC or disk diffusion values, against epidemiological cut-off (ECOFF)
  values. These will be reported in a new column `ecoff`, of class 'sir'
  and coded as 'R' (nonwildtype) or 'S' (wildtype).

- rename_cols:

  A logical value (default is TRUE). If `TRUE`, the function will rename
  the provided columns (specified by `ab_col`, `mic_col`, `disk_col`,
  `species_col`, `id_col`) to the default names expected by AMRgen
  functions ('drug_agent', 'mic', 'disk', 'spp_pheno', 'id'), to match
  those output by the other
  [`import_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ast.md)
  functions.

## Value

A data frame with the processed AST data, including additional columns:

## Examples

``` r
if (FALSE) { # \dontrun{
# import and process AST data from EBI, write formatted data to file for later use
pheno <- import_ebi_ast("EBI_AMR_data.csv.gz")
write_tsv(pheno, file="EBI_AMR_data_processed.tsv.gz", 
            interpret_eucast = TRUE, interpret_ecoff = TRUE)

# read stored data and format the columns to the correct classes
pheno <- format_ast("EBI_AMR_data_processed.tsv.gz")

# read in unprocessed E. coli AST data from non-standard format and interpret
pheno <- format_ast("AMR_data.tsv", sample_col="STRAIN", species="E. coli",
            ab_col="Antibiotic", mic_col="MIC (mg/L)", 
            interpret_eucast=TRUE, interpret_ecoff=TRUE)
} # }
```
