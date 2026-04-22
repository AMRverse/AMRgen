# Summarise a Genotype Table

`summarise_geno()` computes summary information for a genotype table.

## Usage

``` r
summarise_geno(
  geno_table,
  sample_col = "id",
  marker_col = "marker",
  drug_col = "drug",
  class_col = "drug_class",
  gene_col = "gene",
  variation_col = "variation type",
  force_ab = FALSE
)
```

## Arguments

- geno_table:

  A tibble or data frame containing genotype data, in the format output
  by [import_amrfp](https://amrgen.org/reference/import_amrfp.md).

- sample_col:

  Character. Name of the column containing sample identifiers. Default
  is `"id"`.

- marker_col:

  Character. Name of the column containing marker identifiers. Default
  is `"marker"`.

- drug_col:

  Character. Name of the column containing drug agent identifiers.
  Default is `"drug"`. If this is of class 'ab' the entries will be
  annotated with their full antibiotic names, converted using
  [AMR::as.ab](https://amr-for-r.org/reference/as.ab.html). If this is
  desired behaviour but the class is not 'ab', set `force_ab=TRUE`.

- class_col:

  Character. Name of the column containing drug class identifiers.
  Default is `"drug_class"`.

- gene_col:

  Character. Name of the column containing gene identifiers. Default is
  `"gene"`.

- variation_col:

  Character. Name of the column containing variation type identifiers.
  Default is `"variation type"`.

- force_ab:

  Logical. If `TRUE`, attempts to convert entries in `drug_col` to
  antibiotic names using
  [AMR::as.ab](https://amr-for-r.org/reference/as.ab.html) even if this
  column is not of class `"ab"` Default is `FALSE`.

## Value

A named list with the following elements:

- `uniques`: A tibble of the number of unique samples, markers, genes,
  drugs, classes and variation types detected in `geno_table`.

- `per_type`: A tibble of unique counts of samples, markers, genes,
  drugs, and classes per variation type.

- `drugs`: A tibble listing the drugs and/or drug classes represented in
  the table, and the associated number of unique markers, unique
  samples, and total hits for each drug/class.

- `markers`: A tibble listing the markers represented in the table, and
  the associated drugs/classes and variation types (if present). Number
  indicates the count of hits detected per marker.

## Details

The function automatically adapts to the presence or absence of columns
in `geno_table`. The `force_ab` parameter allows the addition of full
antibiotic names using the `ab_name()` function even when the first
column is not recognized as an `"ab"` object.

## Examples

``` r
summarise_geno(staph_geno_ebi)
#> $uniques
#> # A tibble: 5 × 2
#>   column     n_unique
#>   <chr>         <int>
#> 1 id             7547
#> 2 marker           22
#> 3 drug             13
#> 4 drug_class        6
#> 5 gene             22
#> 
#> $per_type
#> NULL
#> 
#> $drugs
#> # A tibble: 15 × 6
#>    drug drug_name     drug_class      markers samples  hits
#>    <ab> <chr>         <chr>             <int>   <int> <int>
#>  1 AMK  Amikacin      Aminoglycosides       6    1941  7745
#>  2 APR  Apramycin     Aminoglycosides       1      12    12
#>  3 CLI  Clindamycin   Lincosamides          1       1     8
#>  4 DOX  Doxycycline   Tetracyclines         1       1     8
#>  5 ERY  Erythromycin  Macrolides            1       1     8
#>  6 GEN  Gentamicin    Aminoglycosides       5     981  4321
#>  7 KAN  Kanamycin     Aminoglycosides       6    1963  7861
#>  8 MNO  Minocycline   Tetracyclines         1       1     8
#>  9 SPT  Spectinomycin Other                 2    1475  1671
#> 10 STR1 Streptomycin  Aminoglycosides       3     497   582
#> 11 TGC  Tigecycline   Tetracyclines         2    7538  7576
#> 12 TOB  Tobramycin    Aminoglycosides       5    1004  4445
#> 13 NA   NA            Aminoglycosides       2      43    43
#> 14 NA   NA            Pleuromutilins        1       1     8
#> 15 NA   NA            Tetracyclines         6    7525 11649
#> 
#> $markers
#> # A tibble: 43 × 5
#>    marker                 drug drug_name     drug_class          n
#>    <chr>                  <ab> <chr>         <chr>           <int>
#>  1 aac(6')-Ie             AMK  Amikacin      Aminoglycosides     9
#>  2 aac(6')-Ie             KAN  Kanamycin     Aminoglycosides     9
#>  3 aac(6')-Ie             TOB  Tobramycin    Aminoglycosides     9
#>  4 aac(6')-Ie/aph(2'')-Ia AMK  Amikacin      Aminoglycosides  4272
#>  5 aac(6')-Ie/aph(2'')-Ia GEN  Gentamicin    Aminoglycosides  4272
#>  6 aac(6')-Ie/aph(2'')-Ia KAN  Kanamycin     Aminoglycosides  4272
#>  7 aac(6')-Ie/aph(2'')-Ia TOB  Tobramycin    Aminoglycosides  4272
#>  8 aadD1                  KAN  Kanamycin     Aminoglycosides   124
#>  9 aadD1                  TOB  Tobramycin    Aminoglycosides   124
#> 10 ant(3'')-IIa           SPT  Spectinomycin Other               2
#> # ℹ 33 more rows
#> 
```
