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
#> # A tibble: 3 × 2
#>   column     n_unique
#>   <chr>         <int>
#> 1 marker           22
#> 2 drug_class        6
#> 3 gene             22
#> 
#> $per_type
#> NULL
#> 
#> $drugs
#> # A tibble: 6 × 3
#>   drug_class              markers     n
#>   <chr>                     <int> <int>
#> 1 Aminoglycosides              14 25009
#> 2 Macrolides/lincosamides       1    16
#> 3 Other antibacterials          2  1671
#> 4 Pleuromutilins                1     8
#> 5 Tetracyclines                 8 19225
#> 6 NA                            1    16
#> 
#> $markers
#> # A tibble: 27 × 3
#>    marker                 drug_class               n
#>    <chr>                  <chr>                <int>
#>  1 aac(6')-Ie             Aminoglycosides         27
#>  2 aac(6')-Ie/aph(2'')-Ia Aminoglycosides      17088
#>  3 aadD1                  Aminoglycosides        248
#>  4 ant(3'')-IIa           Aminoglycosides          2
#>  5 ant(3'')-IIa           Other antibacterials     2
#>  6 ant(6)-Ia              Aminoglycosides        442
#>  7 ant(9)-Ia              Other antibacterials  1669
#>  8 aph(2'')-I             Aminoglycosides         16
#>  9 aph(2'')-Ia            Aminoglycosides        144
#> 10 aph(3')-IIIa           Aminoglycosides       6832
#> # ℹ 17 more rows
#> 
```
