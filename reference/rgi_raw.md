# Example Resistance Gene Identifier (RGI) v6.0.6 Genotype Data

Raw RGI v6.0.6 results file (run with `--include_loose`) for 12 genomes
of multiple species, one AMR determinant per row. Includes multiple
species to cover all four model types currently detected by RGI (protein
homolog model, protein variant model, protein overexpression model, and
rRNA gene variant model) Includes Perfect, Strict, and Loose hits to
test `exclude_loose` parameter

## Usage

``` r
rgi_raw
```

## Format

`rgi_raw` A data frame with 21,203 rows and 28 columns:

- `ORF_ID`: Sample identifier

- ...: RGI results columns

## Source

Four colistin-resistant isolates from Bioproject PRJNA966919
<https://www.ncbi.nlm.nih.gov/datasets/genome/?bioproject=PRJNA966919>.

One genome with rRNA gene variant models detected (GCF_000249055.1).
