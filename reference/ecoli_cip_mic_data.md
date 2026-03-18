# EUCAST Reference distribution for Ciprofloxacin in E. coli

Data frame containing EUCAST reference distribution for ciprofloxacin in
E. coli, downloaded using
[get_eucast_mic_distribution](https://amrgen.org/reference/get_eucast_amr_distribution.md).

## Usage

``` r
ecoli_cip_mic_data
```

## Format

A data frame with 19 rows and 2 columns. It provides MIC distributions
from EUCAST in the form of counts per value.

Columns include:

- `mic`: MIC value.

- `count`: Count of samples with this MIC value, downloaded from EUCAST
  (Feb 2026).

## Source

<https://mic.eucast.org/>
