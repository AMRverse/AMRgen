# Import and Process AST Data from files downloaded from the EBI AMR portal website

This function imports an antibiotic susceptibility testing (AST) dataset
that has been downloaded from the EBI AMR portal website
(https://www.ebi.ac.uk/amr/data/?view=experiments) Data downloaded from
the EBI AMR Portal FTP site
(ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/), either
directly or via the function
[`download_ebi()`](https://AMRverse.github.io/AMRgen/reference/download_ebi.md),
is formatted differently and can instead be processed using the
[`import_ebi_ast_ftp()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast_ftp.md)
function.

## Usage

``` r
import_ebi_ast(
  input,
  sample_col = "phenotype-BioSample_ID",
  source = NULL,
  species = NULL,
  ab = NULL,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE
)
```

## Arguments

- input:

  A string representing a dataframe, or a path to an input file,
  containing the AST data in EBI antibiogram format. These files can be
  downloaded from the EBI AMR browser, e.g.
  https://www.ebi.ac.uk/amr/data/?view=experiments

- sample_col:

  A string indicating the name of the column with sample identifiers. If
  `NULL`, assume this is 'phenotype-BioSample_ID'.

- source:

  (optional) A single value to record as the source of these data
  points, e.g. "EBI_browser". By default, the field
  'phenotype-AMR_associated_publications' will be used to indicate the
  source for each row in the input file, but if this is missing or you
  want to override it with a single value for all samples, you may
  provide a source name via this parameter.

- species:

  (optional) Name of the species to use for phenotype interpretation. By
  default, the field 'phenotype-organism' will be assumed to specify the
  species for each row in the input file, but if this is missing or you
  want to override it in the interpretation step, you may provide a
  single species name via this parameter.

- ab:

  (optional) Name of the antibiotic to use for phenotype interpretation.
  By default, the field 'phenotype-antibiotic_name' will be assumed to
  specify the antibiotic for each row in the input file, but if this is
  missing or you want to override it in the interpretation step, you may
  provide a single antibiotic name via this parameter.

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

## Value

A data frame with the processed AST data, including additional columns:

- `id`: The biological sample identifier (renamed from
  `phenotype-BioSample_ID` or specified column).

- `spp_pheno`: The species phenotype, formatted using the `as.mo`
  function.

- `drug_agent`: The antibiotic used in the test, formatted using the
  `as.ab` function.

- `mic`: The minimum inhibitory concentration (MIC) value, formatted
  using the `as.mic` function.

- `disk`: The disk diffusion measurement (in mm), formatted using the
  `as.disk` function.

- `method`: The AST method (e.g., "MIC", "disk diffusion", "Etest",
  "agar dilution").

- `platform`: The AST platform/instrument (e.g., "Vitek", "Phoenix",
  "Sensititre").

- `guideline`: The AST standard recorded in the input file as being used
  for the AST assay.

- `pheno_eucast`: The phenotype newly interpreted against EUCAST human
  breakpoint standards (as S/I/R), based on the MIC or disk diffusion
  data.

- `pheno_clsi`: The phenotype newly interpreted against CLSI human
  breakpoint standards (as S/I/R), based on the MIC or disk diffusion
  data.

- `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R),
  based on the MIC or disk diffusion data.

- `pheno_provided`: The original phenotype interpretation provided in
  the input file.

- `source`: The source of each data point (renamed from the publications
  field in the input file, or replaced with a single value passed in as
  the 'source' parameter).

## Details

This function will process the data, and optionally interpret the
results based on MIC or disk diffusion data. It assumes that the input
file is a tab-delimited text file (e.g., TSV) or CSV (which may be
compressed) and parses relevant columns (antibiotic names, species
names, MIC or disk data) into suitable classes using the AMR package. It
optionally can use the AMR package to interpret susceptibility phenotype
(SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or
ECOFF).

## Examples

``` r
if (FALSE) { # \dontrun{
# import without re-interpreting resistance
pheno <- import_ebi_ast("EBI_AMR_data.csv.gz")
head(pheno)

# import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
pheno <- import_ebi_ast("EBI_AMR_data.csv.gz", interpret_eucast = TRUE, interpret_ecoff = TRUE)
head(pheno)
} # }
```
