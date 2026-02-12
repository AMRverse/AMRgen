# Import and process antimicrobial phenotype data retrieved from NCBI BioSamples

This function will import antibiotic susceptibility testing (AST) data
suitable for downstream use with AMRgen analysis functions. The expected
input is phenotype data retrieved from NCBI BioSample database via the
function
[`download_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/download_ncbi_ast.md).
Note that files downloaded from the NCBI AST web browser
<https://www.ncbi.nlm.nih.gov/pathogens/ast> are formatted differently
and can be imported with the function
[`import_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ncbi_ast.md).

## Usage

``` r
import_ncbi_biosample(
  input,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE
)
```

## Arguments

- input:

  A string representing the input dataframe, or a path to an input file,
  to be processed.

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

## Examples

``` r
if (FALSE) { # \dontrun{
# Download Klebsiella quasipneumoniae data, filter to amikacin
ast <- download_ncbi_ast(
  "Klebsiella quasipneumoniae",
  antibiotic = "amikacin"
)

# Reformat to simplify use with AMRgen functions
ast <- import_ncbi_biosample(ast, interpret_eucast = T)
} # }
```
