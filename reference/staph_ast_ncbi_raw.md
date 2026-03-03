# S. aureus Example of Raw Phenotype Data Downloaded from NCBI BioSamples via Entrez API

Phenotypes sourced from NCBI Biosamples using the
[download_ncbi_ast](https://AMRverse.github.io/AMRgen/reference/download_ncbi_ast.md)
function without reformating.

## Usage

``` r
staph_ast_ncbi_raw
```

## Format

`staph_ast_ncbi_raw` A data frame with 143 rows and 13 columns
representing all Staphylococcus aureus phenotyping results for amikacin
and doxycycline.

Columns include:

- `id`: Sample identifier.

- `Antibiotic`: Antibiotic name.

- `Resistance phenotype`: S/I/R phenotypes as downloaded from NCBI.

- ...: Additional data columns from NCBI.

## Source

<https://www.ncbi.nlm.nih.gov/pathogens/ast>
