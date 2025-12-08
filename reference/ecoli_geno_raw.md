# E. coli Genotype Example Data

Genotypes called using AMRFinderPlus (v3.12.8, DB 2024-01-31.1), sourced
from the AllTheBacteria project.

## Usage

``` r
ecoli_geno_raw
```

## Format

### `ecoli_geno_raw` A data frame with 45228 rows and 28 columns representing genotyping results from AMRFinderPlus.

Columns include:

- `Name`: Sample identifier.

- `Gene symbol`: Gene symbol in NCBI RefGene.

- `Hierarchy node`: Node in NCBI hierarchy.

- `Class`, `Subclass`: Drug class(es) associated with the marker (from
  NCBI RefGene).

- `% Coverage of reference sequence`,
  `% Identity to reference sequence`, `Accession of closest sequence`:
  Sequence match information.

- ...: Additional metadata columns from the AMRFinderPlus output.

## Source

<https://github.com/ncbi/amr/wiki/Interpreting-results>
