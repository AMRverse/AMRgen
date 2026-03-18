# Import and Process Resistance Gene Identifier (RGI) Results

This function imports and processes genotyping results from the
Resistance Gene Identifier (RGI, https://github.com/arpcard/rgi),
extracting antimicrobial resistance determinants and mapping them to
standardised drug classes/antibiotics.

## Usage

``` r
import_rgi(
  input_table,
  orf_id_col = "ORF_ID",
  sample_id_sep = ".fasta.txt:",
  model_col = "Model_type",
  antibiotic_col = "Antibiotic",
  class_col = "Drug Class",
  exclude_loose = TRUE,
  rgi_short_name = rgi_short_name_table,
  rgi_drugs = rgi_drugs_table
)
```

## Arguments

- input_table:

  A character string specifying a dataframe or path to the RGI results
  table (TSV format).

- orf_id_col:

  A character string specifying the column that identifies open reading
  frame ID (ORF_ID) in the dataset (default `ORF_ID`). This column
  includes the sample ID and the contig / genomic location and is a
  default output of RGI.

- sample_id_sep:

  A character string specifying the separator by which the sample ID is
  separated from the remaining text in `ORF_ID` (Default: `.fasta.txt:`)
  . For example: in the `ORF_ID` column, "SAMEA3498968.fasta.txt:1_96 \#
  109511 \# 110635....", the sample ID separator is `.fasta.txt:`.

- model_col:

  A character string specifying the column that identifies model type
  identified by RGI (default `Model_type`).

- antibiotic_col:

  Character string specifying the antibiotic column (default
  `Antibiotic`).

- class_col:

  Character string specifying the drug class column (default
  `Drug Class`).

- exclude_loose:

  Logical indicating whether to exclude Loose hits (AMR markers that
  fall below a curated bitscore cutoff as defined by CARD/RGI). Default
  `TRUE`, which excludes Loose hits.

- rgi_short_name:

  A tibble containing a reference table mapping model IDs (from
  CARD/RGI) to shortened model names as provided by CARD
  (https://card.mcmaster.ca/download in aro_index.tsv). Defaults to
  `rgi_short_name_table`, which is provided internally.

- rgi_drugs:

  A tibble containing a reference table mapping CARD drug class / drug
  agents to standardised drug classes/names. Defaults to
  `rgi_drugs_table`, which is provided internally.

## Value

A tibble containing the processed AMR determinants and drug classes that
is AMRgen compatible. The output retains the original columns from the
RGI output along with the newly mapped variables.

## Details

The function performs the following steps:

- Reads the RGI output table.

- Transforms RGI output into long form (i.e., one AMR determinant AND
  drug class / antibiotic per row).

- Maps CARD drug classes and antibiotics to standardised names. This
  processing ensures compatibility with downstream AMRgen analysis
  workflows.

## Examples

``` r
# example RGI data (including Perfect, Strict, and Loose hits)
rgi_raw
#> # A tibble: 21,203 × 29
#>    ORF_ID Contig Start  Stop Orientation Cut_Off Pass_Bitscore Best_Hit_Bitscore
#>    <chr>  <chr>  <dbl> <dbl> <chr>       <chr>           <dbl>             <dbl>
#>  1 GCA_0… JASER…     2   283 +           Loose            1150              24.3
#>  2 GCA_0… JASER…   367  1662 +           Loose             700             172. 
#>  3 GCA_0… JASER…  1666  1989 -           Loose             500              23.5
#>  4 GCA_0… JASER…  2031  3386 -           Loose            1900              29.3
#>  5 GCA_0… JASER…  3507  6158 -           Loose             275              30.4
#>  6 GCA_0… JASER…  6961  7386 -           Loose             910              25.4
#>  7 GCA_0… JASER…  7590  8675 +           Loose             450              45.8
#>  8 GCA_0… JASER…  9735 10118 +           Loose             600              25.8
#>  9 GCA_0… JASER… 10164 11495 -           Loose             500              26.6
#> 10 GCA_0… JASER… 11627 12364 +           Loose             400              37.7
#> # ℹ 21,193 more rows
#> # ℹ 21 more variables: Best_Hit_ARO <chr>, Best_Identities <dbl>, ARO <dbl>,
#> #   Model_type <chr>, SNPs_in_Best_Hit_ARO <chr>, Other_SNPs <chr>,
#> #   `Drug Class` <chr>, `Resistance Mechanism` <chr>, `AMR Gene Family` <chr>,
#> #   Predicted_DNA <chr>, Predicted_Protein <chr>, CARD_Protein_Sequence <chr>,
#> #   `Percentage Length of Reference Sequence` <dbl>, ID <chr>, Model_ID <dbl>,
#> #   Nudged <lgl>, Note <lgl>, Hit_Start <dbl>, Hit_End <dbl>, …

# import using sample_id_sep=`_genomic.fna.txt:` and include Loose hits
rgi <- import_rgi(rgi_raw, sample_id_sep = "_genomic.fna.txt:", exclude_loose = FALSE)

# example RGI data from EuSCAPE project (including only Perfect and Strict hits)
rgi_EuSCAPE_raw
#> # A tibble: 59,402 × 26
#>    ORF_ID                 Contig  Start   Stop Orientation Cut_Off Pass_Bitscore
#>    <chr>                   <dbl>  <dbl>  <dbl> <chr>       <chr>           <dbl>
#>  1 SAMEA3498968.fasta.tx…      1 109511 110635 -           Strict            700
#>  2 SAMEA3498968.fasta.tx…      1 237710 238072 +           Perfect           150
#>  3 SAMEA3498968.fasta.tx…      1 238059 238388 +           Perfect           150
#>  4 SAMEA3498968.fasta.tx…      1 278325 279185 -           Perfect           550
#>  5 SAMEA3498968.fasta.tx…      1 299277 299651 -           Strict            230
#>  6 SAMEA3498968.fasta.tx…      2 120780 123893 -           Strict           1900
#>  7 SAMEA3498968.fasta.tx…      2 379306 380745 +           Perfect           900
#>  8 SAMEA3498968.fasta.tx…      2 437973 441050 -           Strict           1800
#>  9 SAMEA3498968.fasta.tx…      2 441051 444173 -           Strict           1800
#> 10 SAMEA3498968.fasta.tx…      3 347739 348971 +           Strict            700
#> # ℹ 59,392 more rows
#> # ℹ 19 more variables: Best_Hit_Bitscore <dbl>, Best_Hit_ARO <chr>,
#> #   Best_Identities <dbl>, ARO <dbl>, Model_type <chr>,
#> #   SNPs_in_Best_Hit_ARO <chr>, Other_SNPs <chr>, `Drug Class` <chr>,
#> #   `Resistance Mechanism` <chr>, `AMR Gene Family` <chr>,
#> #   `Percentage Length of Reference Sequence` <dbl>, ID <chr>, Model_ID <dbl>,
#> #   Nudged <lgl>, Note <lgl>, Hit_Start <dbl>, Hit_End <dbl>, …

# import using defaults (sample_id_sep=`.fasta.txt:`, exclude_loose = `TRUE`)
import_rgi(rgi_EuSCAPE_raw)
#> # A tibble: 292,972 × 33
#>    id        marker mutation drug_agent drug_class `variation type` marker.label
#>    <chr>     <chr>  <chr>    <chr>      <chr>      <chr>            <chr>       
#>  1 SAMEA349… Klebs… NA       FOX        Cephalosp… Gene presence d… Kpne_OmpK37…
#>  2 SAMEA349… Klebs… NA       CTX        Cephalosp… Gene presence d… Kpne_OmpK37…
#>  3 SAMEA349… Klebs… NA       ERY        Macrolides Gene presence d… Kpne_KpnE   
#>  4 SAMEA349… Klebs… NA       STR1       Aminoglyc… Gene presence d… Kpne_KpnE   
#>  5 SAMEA349… Klebs… NA       TCY        Tetracycl… Gene presence d… Kpne_KpnE   
#>  6 SAMEA349… Klebs… NA       FEP        Cephalosp… Gene presence d… Kpne_KpnE   
#>  7 SAMEA349… Klebs… NA       CRO        Cephalosp… Gene presence d… Kpne_KpnE   
#>  8 SAMEA349… Klebs… NA       RIF        Rifamycins Gene presence d… Kpne_KpnE   
#>  9 SAMEA349… Klebs… NA       COL        Polymyxins Gene presence d… Kpne_KpnE   
#> 10 SAMEA349… Klebs… NA       COL        Polymyxins Gene presence d… Kpne_KpnE   
#> # ℹ 292,962 more rows
#> # ℹ 26 more variables: ORF_ID <chr>, Contig <dbl>, Start <dbl>, Stop <dbl>,
#> #   Orientation <chr>, Cut_Off <chr>, Pass_Bitscore <dbl>,
#> #   Best_Hit_Bitscore <dbl>, Best_Hit_ARO <chr>, Best_Identities <dbl>,
#> #   ARO <dbl>, Model_type <chr>, Other_SNPs <chr>, `Drug Class` <chr>,
#> #   `Resistance Mechanism` <chr>, `AMR Gene Family` <chr>,
#> #   `Percentage Length of Reference Sequence` <dbl>, ID <chr>, …
```
