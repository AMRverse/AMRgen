# Import and Process AST Data files retrieved from the EBI AMR portal FTP site

This function will import antibiotic susceptibility testing (AST) data
suitable for downstream use with AMRgen analysis functions. The expected
input is phenotype data retrieved from the EBI AMR Portal FTP site
(ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) either directly
or via the function
[`download_ebi()`](https://AMRverse.github.io/AMRgen/reference/download_ebi.md).
Note that files downloaded from the EBI AMR Portal web browser
(https://www.ebi.ac.uk/amr/data/?view=experiments) are formatted
differently and can be imported with the function
[`import_ebi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast.md)

## Usage

``` r
import_ebi_ast_ftp(
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
# download Salmonella phenotype data from EBI
pheno_salmonella <- download_ebi(genus="Salmonella")

# reformat to simplify use with AMRgen functions
pheno_salmonella <- import_ebi_ast_ftp(pheno_salmonella)
} # }
```
