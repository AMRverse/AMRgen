# Example Salmonella Genotype-Phenotype Data

Raw genotype-phenotype data for *Salmonella enterica* genomes, one row
per sample.

## Usage

``` r
salm_raw
```

## Format

`salm_raw` A data frame with 115 rows and 7 columns:

- `Sample`: Sample identifier

- `Source`, `Serovar`: Non-AMR related information about each isolate

- `CpL_Genotype`: List of genotypic parkers (separated by `;`)

- `Ciprofloxacin`, `Levofloxacin`, `Moxifloxacin`: MIC value for each
  drug
