# S. aureus Example of Imported NCBI Phenotype Data

Phenotypes sourced from NCBI Biosamples using the
[download_ncbi_ast](https://amrgen.org/reference/download_ncbi_ast.md)
function and imported to AMRgen phenotype table format.

## Usage

``` r
staph_ast_ncbi
```

## Format

`staph_ast_ncbi` A data frame with 143 rows and 19 columns representing
all Staphylococcus aureus phenotyping results for amikacin and
doxycycline downloaded from NCBI using
[download_ncbi_ast](https://amrgen.org/reference/download_ncbi_ast.md),
imported into AMRgen format using
[import_ast](https://amrgen.org/reference/import_ast.md).

Columns include:

- `id`: Sample identifier.

- `drug_agent`: Antibiotic identifier, as class 'ab'.

- `mic`: MIC data, as class 'mic'.

- `disk`: Disk diffusion zone diameter data, as class 'disk'.

- `pheno_provided`, `pheno_eucast`: S/I/R phenotypes as downloaded from
  NCBI, and as re-interpreted from mic/disk measures against EUCAST 2024
  breakpoints.

- ...: Additional data columns from NCBI.

## Source

<https://www.ncbi.nlm.nih.gov/pathogens/ast>
