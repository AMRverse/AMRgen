# E. coli NCBI AST Example Data

A subset of E. coli phenotype data from the NCBI AST browser.

## Usage

``` r
ecoli_ast_raw
```

## Format

### `ecoli_ast_raw` A data frame with 10 rows and 21 columns representing unprocessed data from the NCBI AST browser.

Columns include:

- `#BioSample`: Sample identifier.

- `Scientific name`: Species identifier.

- `Antibiotic`: Antibiotic name.

- `Testing standard`: Interpretation standard (EUCAST or CLSI).

- `Measurement sign`: Measurement sign (\>, \<, =, etc.) relating to MIC
  measurement.

- `MIC (mg/L)`: Minimum inhibitory concentration.

- `Disk diffusion (mm)`: Disk diffusion zone.

- `Resistance phenotype`: Resistance call (SIR) as submitted.

- ...: Additional metadata columns from the NCBI AST export.

## Source

<https://www.ncbi.nlm.nih.gov/pathogens/ast>
