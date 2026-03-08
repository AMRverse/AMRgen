# Summarise a Genotype Table

`summarise_geno()` computes summary information for a genotype table.

## Usage

``` r
summarise_geno(
  geno_table,
  sample_col = "Name",
  marker_col = "marker",
  drug_col = "drug_agent",
  class_col = "drug_class",
  gene_col = "gene",
  variation_col = "variation type",
  force_ab = FALSE
)
```

## Arguments

- geno_table:

  A tibble or data frame containing genotype data, in the format output
  by
  [import_amrfp](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md).

- sample_col:

  Character. Name of the column containing sample identifiers. Default
  is `"Name"`.

- marker_col:

  Character. Name of the column containing marker identifiers. Default
  is `"marker"`.

- drug_col:

  Character. Name of the column containing drug agent identifiers.
  Default is `"drug_agent"`. If this is of class 'ab' the entries will
  be annotated with their full antibiotic names, converted using
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

- uniques:

  A tibble of the number of unique samples, markers, genes, drugs,
  classes and variation types detected in `geno_table`.

- pertype:

  A tibble of unique counts of samples, markers, genes, drugs, and
  classes per variation type.

- drugs:

  A tibble listing the drugs and/or drug classes represented in the
  table, and the associated number of unique markers, unique samples,
  and total hits for each drug/class.

- markers:

  A tibble listing the markers represented in the table, and the
  associated drugs/classes and variation types (if present). Number
  indicates the count of hits detected per marker.

## Details

The function automatically adapts to the presence or absence of columns
in `geno_table`. The `force_ab` parameter allows the addition of full
antibiotic names using the `ab_name()` function even when the first
column is not recognized as an `"ab"` object.

## Examples

``` r
geno_table <- import_amrfp(ecoli_geno_raw)
summarise_geno(geno_table)
#> $uniques
#> # A tibble: 6 × 2
#>   column         n_unique
#>   <chr>             <int>
#> 1 Name               5258
#> 2 marker              244
#> 3 drug_agent           35
#> 4 drug_class           26
#> 5 gene                196
#> 6 variation type        5
#> 
#> $pertype
#> # A tibble: 5 × 6
#>   `variation type`                Name marker drug_agent drug_class  gene
#>   <chr>                          <int>  <int>      <int>      <int> <int>
#> 1 Gene presence detected          5258    164         22         17   164
#> 2 Inactivating mutation detected   615     42         15         14    42
#> 3 Nucleotide variant detected       57      2          3          3     1
#> 4 Promoter variant detected         93      4          1          1     1
#> 5 Protein variant detected        4920     65         18         16    21
#> 
#> $drugs
#> # A tibble: 44 × 6
#>    drug_agent antibiotic                  drug_class       markers samples  hits
#>    <ab>       <chr>                       <chr>              <int>   <int> <int>
#>  1 AMC        Amoxicillin/clavulanic acid Aminopenicillins       2      57    57
#>  2 AMK        Amikacin                    Aminoglycosides        6     176   180
#>  3 AMP        Ampicillin                  Aminopenicillins       6     749   749
#>  4 APR        Apramycin                   Aminoglycosides        1      98    98
#>  5 ATM        Aztreonam                   Monobactams            2      39    39
#>  6 AZM        Azithromycin                Macrolides             4     472   478
#>  7 BLM        Bleomycin                   Glycopeptides          2      40    40
#>  8 CHL        Chloramphenicol             Phenicols             15    1121  1181
#>  9 CLI        Clindamycin                 Lincosamides           1      26    26
#> 10 CLR        Clarithromycin              Macrolides             1       1     1
#> # ℹ 34 more rows
#> 
#> $markers
#> # A tibble: 349 × 6
#>    marker      drug_agent antibiotic  drug_class      `variation type`         n
#>    <chr>       <ab>       <chr>       <chr>           <chr>                <int>
#>  1 aac(2')-IIa KAS        Kasugamycin Aminoglycosides Gene presence detec…     1
#>  2 aac(3)-II   GEN        Gentamicin  Aminoglycosides Inactivating mutati…     1
#>  3 aac(3)-IId  GEN        Gentamicin  Aminoglycosides Gene presence detec…   204
#>  4 aac(3)-IId  GEN        Gentamicin  Aminoglycosides Inactivating mutati…     1
#>  5 aac(3)-IIe  GEN        Gentamicin  Aminoglycosides Gene presence detec…   120
#>  6 aac(3)-IIg  GEN        Gentamicin  Aminoglycosides Gene presence detec…     2
#>  7 aac(3)-IVa  APR        Apramycin   Aminoglycosides Gene presence detec…    98
#>  8 aac(3)-IVa  GEN        Gentamicin  Aminoglycosides Gene presence detec…    98
#>  9 aac(3)-IVa  TOB        Tobramycin  Aminoglycosides Gene presence detec…    98
#> 10 aac(3)-Ib   GEN        Gentamicin  Aminoglycosides Gene presence detec…     1
#> # ℹ 339 more rows
#> 
```
