# Import and Process Kleborate Results

This function imports and processes genotyping results from Kleborate
(https://github.com/klebgenomics/Kleborate), extracting antimicrobial
resistance determinants and mapping them to standardised drug classes.

## Usage

``` r
import_kleborate(
  input_table,
  sample_col = "strain",
  kleborate_class_table = kleborate_classes,
  hgvs = TRUE
)
```

## Arguments

- input_table:

  A character string specifying a dataframe or path to the Kleborate
  results table (TSV format).

- sample_col:

  A character string specifying the column that identifies samples in
  the dataset (default `strain`).

- kleborate_class_table:

  A tibble containing a reference table mapping Kleborate drug class
  column names (`Kleborate_Class`) to standardised drug classes
  (`drug_class`). Defaults to `kleborate_classes`, which is provided
  internally.

- hgvs:

  Logical indicating whether to expect mutations in [HGVS
  nomenclature](https://hgvs-nomenclature.org/stable/) syntax (used in
  Kleborate releases since v3.1.3). Default `TRUE`, which expects
  mutations formatted as e.g. "GyrA:p.S83F". Set to `FALSE` if your
  results were generated using older versions where mutations were
  formatted as e.g. "GyrA_83F".

## Value

A data frame with the processed genotype data, with harmonised gene
names, mapped drug agents, and drug classes which can be used for other
functions of the ARMgen package:

- `id`: The sample identifier (`character`).

- `marker`: The name of the genotype marker as it appears in the input
  (e.g. `GyrA:p.S83F` for recent versions of Kleborate, or `GyrA-83F`
  for earlier versions not using [HGVS
  nomenclature](https://hgvs-nomenclature.org/stable/)) (`character`).

- `gene`: The gene identifier (`character`).

- `mutation`: The mutation detected within the gene, converted to [HGVS
  nomenclature](https://hgvs-nomenclature.org/stable/) syntax (e.g.
  `Ser83Phe`) (`character`).

- `drug_class`: Name of the antibiotic group associated with the
  genotype marker, compatible with AMR pkg, drawn from the Kleborate
  column in which the marker was reported (`character`).

- `drug`: Values are recorded as `NA` as Kleborate doesn't report
  markers assigned to individual drug level.

- `variation type`: Type of variation, e.g. `Gene presence detected`,
  `Protein variant detected`, `Nucleotide variant detected`,
  `Inactivating mutation detected`. ... Other fields specific to the
  input file

## Details

The function performs the following steps:

- Reads the Kleborate output table.

- Transforms Kleborate output into long form (i.e., one AMR determinant
  per row).

- Maps Kleborate drug classes to standardised drug class names
  recognised by the AMR pkg. This processing ensures compatibility with
  downstream AMRgen analysis workflows.

## Examples

``` r
# example Kleborate data from EUSCAPE project
kleborate_raw
#> # A tibble: 1,490 × 122
#>    strain    species species_match contig_count    N50 largest_contig total_size
#>    <chr>     <chr>   <chr>                <dbl>  <dbl>          <dbl>      <dbl>
#>  1 SAMEA349… Klebsi… strong                 141 230759         470757    5578320
#>  2 SAMEA349… Klebsi… strong                  88 370309         938079    5384685
#>  3 SAMEA349… Klebsi… strong                  90 238750         529125    5446454
#>  4 SAMEA349… Klebsi… strong                 144 207582         663698    5574298
#>  5 SAMEA349… Klebsi… strong                 142 263498         678692    5486238
#>  6 SAMEA349… Klebsi… strong                  79 285199         991412    5529803
#>  7 SAMEA349… Klebsi… strong                 280 178980         585359    5817055
#>  8 SAMEA349… Klebsi… strong                 108 209418         517450    5379124
#>  9 SAMEA349… Klebsi… strong                 134 371444         984005    5558705
#> 10 SAMEA349… Klebsi… strong                 142 197944         636773    5497421
#> # ℹ 1,480 more rows
#> # ℹ 115 more variables: GC_content <dbl>, ambiguous_bases <chr>,
#> #   QC_warnings <chr>, ST <chr>, gapA <dbl>, infB <dbl>, mdh <dbl>, pgi <dbl>,
#> #   phoE <dbl>, rpoB <dbl>, tonB <dbl>, YbST <chr>, Yersiniabactin <chr>,
#> #   ybtS <chr>, ybtX <chr>, ybtQ <chr>, ybtP <chr>, ybtA <chr>, irp2 <chr>,
#> #   irp1 <chr>, ybtU <chr>, ybtT <chr>, ybtE <chr>, fyuA <chr>,
#> #   spurious_ybt_hits <chr>, CbST <chr>, Colibactin <chr>, clbA <chr>, …

# import first few rows of this data frame and parse it to standard genotype table format
kleborate_geno <- import_kleborate(kleborate_raw %>% head(n = 10))

# parse the output of an older version of Kleborate (v3.1.3) before
# HGVS syntax was introduced for mutations
kleborate_geno <- import_kleborate(kleborate_raw_v313 %>% head(n = 10), hgvs = FALSE)
```
