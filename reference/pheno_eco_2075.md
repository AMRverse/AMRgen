# E. coli AST data from Mills et al 2022

Phenotyping data published in Mills et al, Genome Medicine (2022)
14:147, downloaded via EBI AMR portal, and imported to AMRgen phenotype
table format. Corresponding genotype data is available in
`geno_eco_2075`.

## Usage

``` r
pheno_eco_2075
```

## Format

`pheno_eco_2075` A data frame with 37350 rows and 37 columns
representing MIC results for 2075 E. coli isolates tested against 18
drug agents, using BD Phoenix.

Columns include:

- `id`: Sample identifier (BioSample).

- `drug`: Antibiotic identifier, as class 'ab'.

- `mic`: MIC data, as class 'mic'.

- `pheno_provided`: S/I/R phenotypes as downloaded from EBI.

- ...: Additional data columns from EBI.

## Source

<https://www.ebi.ac.uk/amr>
