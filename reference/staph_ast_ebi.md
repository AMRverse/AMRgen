# S. aureus Example of Imported EBI Phenotype Data

Phenotypes sourced from EBI AMR Portal using the
[download_ebi](https://AMRverse.github.io/AMRgen/reference/download_ebi.md)
function and imported to AMRgen phenotype table format.

## Usage

``` r
staph_ast_ebi
```

## Format

### `staph_ast_ebi` A data frame with 218 rows and 46 columns representing phenotyping results from EBI, imported into AMRgen format using [import_ast](https://AMRverse.github.io/AMRgen/reference/import_ast.md).

Columns include:

- `id`: Sample identifier.

- `drug_agent`: Antibiotic identifier, as class 'ab'.

- `mic`: MIC data, as class 'mic'.

- `disk`: Disk diffusion zone diameter data, as class 'disk'.

- `pheno_provided`, `pheno_eucast`, `pheno_clsi`, `ecoff`: S/I/R
  phenotypes as downloaded from EBI, and as re-interpreted from mic/disk
  measures against EUCAST 2024 breakpoints using
  [AMR::as.ab](https://amr-for-r.org/reference/as.ab.html).

- ...: Additional data columns from EBI.

## Source

<https://www.ebi.ac.uk/amr>
