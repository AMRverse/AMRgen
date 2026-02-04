# Interpret AST data in a standard format tibble

This function applies human EUCAST or CLSI breakpoints, and/or ECOFF, to
interpret AST data.

## Usage

``` r
interpret_ast(
  ast,
  interpret_ecoff = TRUE,
  interpret_eucast = TRUE,
  interpret_clsi = TRUE,
  species = NULL,
  ab = NULL
)
```

## Arguments

- ast:

  A tibble containing the AST measures in standard AMRgen format, as
  output by `import_ast`. It must contain assay measurements in columns
  'mic' (class mic) and/or 'disk'. Interpretation requires an organism
  (column 'spp_pheno' of class 'mo', or a single value passed via the
  'species' parameter) and an antibiotic (column 'drug_agent' of class
  'ab', or a single value passed via the 'ab' parameter).

- interpret_ecoff:

  A logical value (default is FALSE). If `TRUE`, the function will
  interpret the wildtype vs nonwildtype status for each row based on the
  MIC or disk diffusion values, against epidemiological cut-off (ECOFF)
  values. These will be reported in a new column `ecoff`, of class 'sir'
  and coded as 'R' (nonwildtype) or 'S' (wildtype).

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

- species:

  (optional) Name of the species to use for phenotype interpretation. By
  default, the spp_pheno field in the input file will be assumed to
  specify the species for each sample, but if this is missing or you
  want to override it in the interpretation step, you may provide a
  single species name via this parameter.

- ab:

  (optional) Name of the antibiotic to use for phenotype interpretation.
  By default, the drug_agent field in the input file will be assumed to
  specify the antibiotic for each sample, but if this is missing or you
  want to override it in the interpretation step, you may provide a
  single antibiotic name via this parameter.

## Value

A copy of the input table, with additional columns:

- `pheno_eucast`: The phenotype newly interpreted against EUCAST human
  breakpoint standards (as S/I/R), based on the MIC or disk diffusion
  data.

- `pheno_clsi`: The phenotype newly interpreted against CLSI human
  breakpoint standards (as S/I/R), based on the MIC or disk diffusion
  data.

- `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R),
  based on the MIC or disk diffusion data.

- `spp_pheno`: The species phenotype, formatted using the `as.mo`
  function (either taken from the input table, or the single value
  specified by 'species' parameter).

- `drug_agent`: The antibiotic used in the test, formatted using the
  `as.ab` function (either taken from the input table, or the single
  value specified by 'ab' parameter)..

## Examples

``` r
# import without re-interpreting resistance
pheno <- import_ncbi_ast(ecoli_ast_raw)
#> Warning: Expected AST method column 'Laboratory typing method' not found in input
#> Warning: Expected column 'BioProject' not found in input
head(pheno)
#> # A tibble: 6 × 29
#>   id           drug_agent     mic  disk guideline method platform pheno_provided
#>   <chr>        <ab>         <mic> <dsk> <chr>     <chr>  <chr>    <sir>         
#> 1 SAMN36015110 CIP        <128.00    NA CLSI      MIC    NA         R           
#> 2 SAMN11638310 CIP         256.00    NA CLSI      MIC    NA         R           
#> 3 SAMN05729964 CIP          64.00    NA CLSI      Etest  Etest      R           
#> 4 SAMN10620111 CIP         >=4.00    NA CLSI      MIC    NA         R           
#> 5 SAMN10620168 CIP         >=4.00    NA CLSI      MIC    NA         R           
#> 6 SAMN10620104 CIP         <=0.25    NA CLSI      MIC    NA         S           
#> # ℹ 21 more variables: spp_pheno <mo>, `Organism group` <chr>,
#> #   `Scientific name` <chr>, `Isolation type` <chr>, Location <chr>,
#> #   `Isolation source` <chr>, Isolate <chr>, Antibiotic <chr>,
#> #   `Resistance phenotype` <chr>, `Measurement sign` <chr>, `MIC (mg/L)` <dbl>,
#> #   `Disk diffusion (mm)` <lgl>, `Laboratory typing platform` <chr>,
#> #   Vendor <chr>, `Laboratory typing method version or reagent` <chr>,
#> #   `Testing standard` <chr>, `Create date` <dttm>, pheno_clsi_mic <sir>, …

# interpret phenotypes
pheno <- interpret_ast(pheno)
#> Warning: There was 1 warning in `mutate()`.
#> ℹ In argument: `across(...)`.
#> Caused by warning:
#> ! Some MICs were converted to the nearest higher log2 level, following the
#> CLSI interpretation guideline.

if (FALSE) { # \dontrun{
pheno <- read_csv("AST.csv") %>% 
  mutate(drug_agent=as.ab(antibiotic)) %>% # convert antibiotic field to 'drug_agent' of class 'ab'
  mutate(mic=paste0(sign,MIC)) %>% 
  mutate(mic=as.mic(mic)) # create a single 'mic' column of class 'mic'
  
pheno <- interpret_ast(pheno, species="Escherichia coli")
} # }
```
