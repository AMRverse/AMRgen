# NCBI AST for Escherichia coli tested against chloramphenicol

NCBI Antibiotic Susceptibility Test (AST) Browser data for Escherichia
coli tested against chloramphenicol.

## Usage

``` r
NCBI_Ecoli_pheno_chl
```

## Format

`NCBI_Ecoli_pheno_chl` A data frame with 6,859 rows and 17 columns:

- `id`: Sample identifier, imported from the `BioSample` column in the
  raw input.

- `drug`: Antibiotic code, interpreted from `Antibiotic` using `as.ab`.

- `mic`: Minimum inhibitory concentration, formatted using `as.mic`.

- `disk`: Disk diffusion zone, formatted using `as.disk`.

- `method`, `platform`, `guideline`: Test method and platform and
  interpretation guideline.

- `pheno_provided`: S/I/R interpretation as provided in the raw input.

- `spp_pheno`: Species identifier, interpreted from `Scientific name`
  using `as.mo`, used to interpret `ecoff` and `pheno` columns.

- ...: Additional data columns from the NCBI AST Browser

## Source

<https://www.ncbi.nlm.nih.gov/pathogens/ast#chloramphenicol%20AND%20Escherichia>
