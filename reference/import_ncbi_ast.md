# Import and Process AST Data from an NCBI File

This function imports an antibiotic susceptibility testing (AST)
dataset, processes the data, and optionally interprets the results based
on MIC or disk diffusion data. It assumes that the input file is a
tab-delimited text file (e.g., TSV) or CSV (which may be compressed) and
parses relevant columns (antibiotic names, species names, MIC or disk
data) into suitable classes using the AMR package. It optionally can use
the AMR package to interpret susceptibility phenotype (SIR) based on
EUCAST or CLSI guidelines (human breakpoints and/or ECOFF). If expected
columns are not found warnings will be given, and interpretation may not
be possible.

## Usage

``` r
import_ncbi_ast(
  input,
  sample_col = "BioSample",
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
  containing the AST data in NCBI antibiogram format. These files can be
  downloaded from NCBI AST browser, e.g.
  https://www.ncbi.nlm.nih.gov/pathogens/ast#Pseudomonas%20aeruginosa

- sample_col:

  A string indicating the name of the column with sample identifiers. If
  `NULL`, assume this is 'BioSample'.

- source:

  (optional) A single value to record as the source of these data
  points, e.g. "NCBI_browser". By default, the field 'BioProject' will
  be used to indicate the source for each row in the input file, but if
  this is missing or you want to override it with a single value for all
  samples, you may provide a source name via this parameter.

- species:

  (optional) Name of the species to use for phenotype interpretation. By
  default, the field 'Scientific name' will be assumed to specify the
  species for each row in the input file, but if this is missing or you
  want to override it in the interpretation step, you may provide a
  single species name via this parameter.

- ab:

  (optional) Name of the antibiotic to use for phenotype interpretation.
  By default, the field 'Antibiotic' will be assumed to specify the
  antibiotic for each row in the input file, but if this is missing or
  you want to override it in the interpretation step, you may provide a
  single antibiotic name via this parameter.

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

- `id`: The biological sample identifier (renamed from `BioSample` or
  specified column).

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

- `source`: The source of each data point (renamed from `BioProject` in
  the input file, or replaced with a single value passed in as the
  'source' parameter).

## Examples

``` r
# small example E. coli AST data from NCBI
ecoli_ast_raw
#> # A tibble: 10 × 21
#>    `#BioSample` `Organism group`    `Scientific name` `Isolation type` Location 
#>    <chr>        <chr>               <chr>             <chr>            <chr>    
#>  1 SAMN36015110 E.coli and Shigella Escherichia coli  clinical         India: A…
#>  2 SAMN11638310 E.coli and Shigella Escherichia coli  clinical         India    
#>  3 SAMN05729964 E.coli and Shigella Escherichia coli  clinical         Brazil: …
#>  4 SAMN10620111 E.coli and Shigella Escherichia coli  clinical         USA: Roc…
#>  5 SAMN10620168 E.coli and Shigella Escherichia coli  clinical         USA: Roc…
#>  6 SAMN10620104 E.coli and Shigella Escherichia coli  clinical         USA: Roc…
#>  7 SAMN10620102 E.coli and Shigella Escherichia coli  clinical         USA: Roc…
#>  8 SAMN10620129 E.coli and Shigella Escherichia coli  clinical         USA: Roc…
#>  9 SAMN10620121 E.coli and Shigella Escherichia coli  clinical         USA: Roc…
#> 10 SAMN10620086 E.coli and Shigella Escherichia coli  clinical         USA: Roc…
#> # ℹ 16 more variables: `Isolation source` <chr>, Isolate <chr>,
#> #   Antibiotic <chr>, `Resistance phenotype` <chr>, `Measurement sign` <chr>,
#> #   `MIC (mg/L)` <dbl>, `Disk diffusion (mm)` <lgl>,
#> #   `Laboratory typing platform` <chr>, Vendor <chr>,
#> #   `Laboratory typing method version or reagent` <chr>,
#> #   `Testing standard` <chr>, `Create date` <dttm>, pheno_clsi_mic <sir>,
#> #   pheno_clsi_disk <sir>, ecoff_mic <sir>, ecoff_disk <sir>

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

# import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
pheno <- import_ncbi_ast(ecoli_ast_raw, interpret_eucast = TRUE, interpret_ecoff = TRUE)
#> Warning: Expected AST method column 'Laboratory typing method' not found in input
#> Warning: Expected column 'BioProject' not found in input
head(pheno)
#> # A tibble: 6 × 33
#>   id       drug_agent     mic  disk pheno_eucast ecoff guideline method platform
#>   <chr>    <ab>         <mic> <dsk> <sir>        <sir> <chr>     <chr>  <chr>   
#> 1 SAMN360… CIP        <128.00    NA   NI           NI  CLSI      MIC    NA      
#> 2 SAMN116… CIP         256.00    NA   R            R   CLSI      MIC    NA      
#> 3 SAMN057… CIP          64.00    NA   R            R   CLSI      Etest  Etest   
#> 4 SAMN106… CIP         >=4.00    NA   R            R   CLSI      MIC    NA      
#> 5 SAMN106… CIP         >=4.00    NA   R            R   CLSI      MIC    NA      
#> 6 SAMN106… CIP         <=0.25    NA   S            NI  CLSI      MIC    NA      
#> # ℹ 24 more variables: pheno_provided <sir>, spp_pheno <mo>,
#> #   `Organism group` <chr>, `Scientific name` <chr>, `Isolation type` <chr>,
#> #   Location <chr>, `Isolation source` <chr>, Isolate <chr>, Antibiotic <chr>,
#> #   `Resistance phenotype` <chr>, `Measurement sign` <chr>, `MIC (mg/L)` <dbl>,
#> #   `Disk diffusion (mm)` <lgl>, `Laboratory typing platform` <chr>,
#> #   Vendor <chr>, `Laboratory typing method version or reagent` <chr>,
#> #   `Testing standard` <chr>, `Create date` <dttm>, pheno_clsi_mic <sir>, …
```
