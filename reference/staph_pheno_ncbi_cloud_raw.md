# S. aureus Example of Raw Phenotype Data Downloaded from NCBI via Google Cloud BigQuery

Phenotypes sourced from NCBI via
[query_ncbi_bq_pheno](https://amrgen.org/reference/query_ncbi_bq_pheno.md)
function, without reformatting.

## Usage

``` r
staph_pheno_ncbi_cloud_raw
```

## Format

`staph_pheno_ncbi_cloud_raw` A data frame with 142 rows and 11 columns
representing all Staphylococcus aureus phenotyping results for amikacin
and doxycycline.

Columns include:

- `BioSample`: Sample identifier.

- `Antibiotic`: Antibiotic name.

- `Resistance phenotype`: S/I/R phenotypes as downloaded from NCBI.

- ...: Additional data columns from NCBI.

## Source

<https://www.ncbi.nlm.nih.gov/pathogens/ast>
