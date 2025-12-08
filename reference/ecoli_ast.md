# E. coli NCBI AST Example Data, Re-interpreted with AMR Package

A subset of E. coli phenotype data from the NCBI AST browser.

## Usage

``` r
ecoli_ast
```

## Format

### `ecoli_ast` A data frame with 4170 rows and 10 columns representing reinterpreted data from the NCBI AST browser.

Columns include:

- `id`: Sample identifier, imported from the `#BioSample` column in the
  raw input.

- `drug_agent`: Antibiotic code, interpreted from `Antibiotic` using
  `as.ab`, used to interpret `ecoff` and `pheno` columns.

- `mic`: Minimum inhibitory concentration, formatted using `as.mic`,
  used to interpret `ecoff` and `pheno` columns.

- `disk`: Disk diffusion zone, formatted using `as.disk`, used to
  interpret `ecoff` and `pheno` columns.

- `pheno_clsi`: S/I/R classification according to CLSI, interpreted
  using `as.sir`.

- `ecoff`: WT/NWT classification, interpreted using `as.sir`.

- `guideline`: Interpretation guidelines used to interpret `ecoff` and
  `pheno` columns.

- `method`: Test method, one of: "Etest', "Microscan', "Phoenix',
  "Phoenix NMIC-203 card', "Sensititer', "Sensititre', "Vitek'.

- `pheno_provided`: ??

- `spp_pheno`: Species identifier, interpreted from `Scientific name`
  using `as.mo`, used to interpret `ecoff` and `pheno` columns.

## Source

<https://www.ncbi.nlm.nih.gov/pathogens/ast>
