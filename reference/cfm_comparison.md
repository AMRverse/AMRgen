# N. gonorrhoeae Euro-GASP cefixime data vs EUCAST reference distribution

Minimum inhibitory concentration (MIC) distributions for cefixime in
*Neisseria gonorrhoeae* isolates. EUCAST reference distribution
(retrieved March 2026 using
[`compare_mic_with_eucast()`](https://AMRverse.github.io/AMRgen/reference/get_eucast_amr_distribution.md))
compared with data from three European Gonococcal Antimicrobial
Surveillance Programme (Euro-GASP) genomic surveys (2013, 2018, 2020),
see Harris *et al.* (2018)
<https://doi.org/10.1016/S1473-3099(18)30225-1>, Sánchez-Busó *et al.*
(2022) <https://doi.org/10.1016/S2666-5247(22)00044-1>, and Golparian
*et al.* (2024) <https://doi.org/10.1016/S2666-5247(23)00370-1>.

## Usage

``` r
cfm_comparison
```

## Format

`cfm_comparison` A data frame with 15 rows and 3 columns:

- `value`: MIC value

- `user`: Number of samples with each MIC value, from Euro-GASP data
  (total n=5,361 isolates)

- `eucast`: Number of samples with each MIC value, from EUCAST reference
  data (total n=35,582 isolates)

## Source

https://mic.eucast.org/
