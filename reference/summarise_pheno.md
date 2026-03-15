# Summarise a Phenotype Table

`summarise_pheno()` computes summary information for a phenotype table.

## Usage

``` r
summarise_pheno(
  pheno_table,
  sample_col = "id",
  drug_col = "drug_agent",
  mic_col = "mic",
  disk_col = "disk",
  spp_col = "spp_pheno",
  pheno_cols = NULL,
  method_cols = c("method", "platform", "guideline", "source"),
  force_ab = FALSE
)
```

## Arguments

- pheno_table:

  A tibble or data frame containing phenotype data, in the format output
  by
  [import_ast](https://AMRverse.github.io/AMRgen/reference/import_ast.md).

- sample_col:

  Character. Name of the column containing sample identifiers. Default
  is `"id"`.

- drug_col:

  Character. Name of the column containing drug agent identifiers.
  Default is `"drug_agent"`. If this is of class 'ab' the entries will
  be annotated with their full antibiotic names, converted using
  [AMR::as.ab](https://amr-for-r.org/reference/as.ab.html). If this is
  desired behaviour but the class is not 'ab', set `force_ab=TRUE`.

- mic_col:

  Character. Name of the column containing MIC measurements Default is
  `"mic"`.

- disk_col:

  Character. Name of the column containing disk diffusion zone
  measurements. Default is `"disk"`.

- spp_col:

  Character. Name of the column containing species names. Default is
  `"spp_pheno"`.

- pheno_cols:

  Vector. Vector giving names of columns containing categorical
  phenotype calls (S/I/R or NWT/WT). Default is any columns beginning
  with `"pheno"` or `"ecoff"`.

- method_cols:

  Vector. Vector giving names of columns containing method or source
  information by which to summarise MIC/disk data. Default is
  `c("method", "platform", "guideline", "source")`.

- force_ab:

  Logical. If `TRUE`, attempts to convert entries in `drug_col` to
  antibiotic names using
  [AMR::as.ab](https://amr-for-r.org/reference/as.ab.html) even if this
  column is not of class `"ab"` Default is `FALSE`.

## Value

A named list with the following elements:

- uniques:

  A tibble of the number of unique samples, drugs, organisms, and
  methods detected in `pheno_table`.

- drugs:

  A tibble listing the drugs included in the table, and the associated
  number of samples with MIC measures, disk measures, neither or both,
  for each drug and species.

- details:

  A tibble listing more details of the methods of assay measurements,
  per drug and species.

- pheno_counts_list:

  A list of tibbles, each corresponding to a unique categorical
  phenotype column in the input, indicating the counts of each
  phenotypic category per drug and species.

## Details

The function automatically adapts to the presence or absence of columns
in `pheno_table`. The `force_ab` parameter allows the addition of full
antibiotic names using the `ab_name()` function even when the first
column is not recognized as an `"ab"` object.

## Examples

``` r
summarise_pheno(staph_ast_ebi)
#> No phenotype column names provided via 'pheno_cols'
#> These are needed to summarise counts of phenotype category calls per drug.
#> Relevant columns detected in your input table are: c('pheno_provided','pheno_eucast','pheno_clsi','pheno_eucast_mic','pheno_eucast_disk','pheno_clsi_mic','pheno_clsi_disk','ecoff','ecoff_mic','ecoff_disk')
#> $uniques
#> # A tibble: 7 × 2
#>   column     n_unique
#>   <chr>         <int>
#> 1 id              190
#> 2 drug_agent        2
#> 3 spp_pheno         1
#> 4 method            3
#> 5 platform          2
#> 6 guideline         2
#> 7 source            2
#> 
#> $drugs
#> # A tibble: 2 × 6
#>   drug_agent antibiotic_name spp_pheno              disk   mic  none
#>   <ab>       <chr>           <chr>                 <int> <int> <int>
#> 1 AMK        Amikacin        Staphylococcus aureus     3    11    70
#> 2 DOX        Doxycycline     Staphylococcus aureus    NA    47    87
#> 
#> $details
#> # A tibble: 6 × 10
#>   drug_agent antibiotic_name spp_pheno    method platform guideline source   mic
#>   <ab>       <chr>           <chr>        <chr>  <chr>    <chr>     <chr>  <int>
#> 1 AMK        Amikacin        Staphylococ… broth… Vitek    CLSI      NA         1
#> 2 AMK        Amikacin        Staphylococ… disk … NA       CLSI      NA        NA
#> 3 AMK        Amikacin        Staphylococ… disk … NA       NA        27150…    NA
#> 4 AMK        Amikacin        Staphylococ… NA     NA       NA        NA        10
#> 5 DOX        Doxycycline     Staphylococ… broth… NA       CLSI      NA        29
#> 6 DOX        Doxycycline     Staphylococ… NA     NA       NA        NA        18
#> # ℹ 2 more variables: disk <int>, none <int>
#> 
#> $pheno_counts_list
#> list()
#> 

summarise_pheno(staph_ast_ebi, pheno_cols = c("pheno_provided", "pheno_clsi", "ecoff"))
#> $uniques
#> # A tibble: 7 × 2
#>   column     n_unique
#>   <chr>         <int>
#> 1 id              190
#> 2 drug_agent        2
#> 3 spp_pheno         1
#> 4 method            3
#> 5 platform          2
#> 6 guideline         2
#> 7 source            2
#> 
#> $drugs
#> # A tibble: 2 × 6
#>   drug_agent antibiotic_name spp_pheno              disk   mic  none
#>   <ab>       <chr>           <chr>                 <int> <int> <int>
#> 1 AMK        Amikacin        Staphylococcus aureus     3    11    70
#> 2 DOX        Doxycycline     Staphylococcus aureus    NA    47    87
#> 
#> $details
#> # A tibble: 6 × 10
#>   drug_agent antibiotic_name spp_pheno    method platform guideline source   mic
#>   <ab>       <chr>           <chr>        <chr>  <chr>    <chr>     <chr>  <int>
#> 1 AMK        Amikacin        Staphylococ… broth… Vitek    CLSI      NA         1
#> 2 AMK        Amikacin        Staphylococ… disk … NA       CLSI      NA        NA
#> 3 AMK        Amikacin        Staphylococ… disk … NA       NA        27150…    NA
#> 4 AMK        Amikacin        Staphylococ… NA     NA       NA        NA        10
#> 5 DOX        Doxycycline     Staphylococ… broth… NA       CLSI      NA        29
#> 6 DOX        Doxycycline     Staphylococ… NA     NA       NA        NA        18
#> # ℹ 2 more variables: disk <int>, none <int>
#> 
#> $pheno_counts_list
#> $pheno_counts_list$pheno_provided
#> # A tibble: 2 × 6
#>   drug_agent antibiotic_name spp_pheno                 S     R     I
#>   <ab>       <chr>           <chr>                 <int> <int> <int>
#> 1 AMK        Amikacin        Staphylococcus aureus    12    72    NA
#> 2 DOX        Doxycycline     Staphylococcus aureus   103    20    11
#> 
#> $pheno_counts_list$pheno_clsi
#> # A tibble: 2 × 7
#>   drug_agent antibiotic_name spp_pheno              `NA`     S     I     R
#>   <ab>       <chr>           <chr>                 <int> <int> <int> <int>
#> 1 AMK        Amikacin        Staphylococcus aureus    84    NA    NA    NA
#> 2 DOX        Doxycycline     Staphylococcus aureus    87    41     5     1
#> 
#> $pheno_counts_list$ecoff
#> # A tibble: 2 × 7
#>   drug_agent antibiotic_name spp_pheno                WT   NWT  `NA`    NI
#>   <ab>       <chr>           <chr>                 <int> <int> <int> <int>
#> 1 AMK        Amikacin        Staphylococcus aureus    12     2    70    NA
#> 2 DOX        Doxycycline     Staphylococcus aureus    NA     9    87    38
#> 
#> 
```
