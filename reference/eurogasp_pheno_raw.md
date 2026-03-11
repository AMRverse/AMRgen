# N. gonorrhoeae Euro-GASP Phenotype Data Use Case 1

Minimum inhibitory concentration (MIC) data for *Neisseria gonorrhoeae*
isolates from three European Gonococcal Antimicrobial Surveillance
Programme (Euro-GASP) genomic surveys (2013, 2018, 2020), from Harris
*et al.* (2018) <https://doi.org/10.1016/S1473-3099(18)30225-1>,
Sánchez-Busó *et al.* (2022)
<https://doi.org/10.1016/S2666-5247(22)00044-1>, and Golparian *et al.*
(2024) <https://doi.org/10.1016/S2666-5247(23)00370-1>.

## Usage

``` r
eurogasp_pheno_raw
```

## Format

`eurogasp_pheno_raw` A data frame with 5,361 rows and 5 columns:

- `id`: Sample identifier (ENA run accession).

- `Azithromycin`: MIC value in mg/L for azithromycin (numeric; n =
  5,055).

- `Ciprofloxacin`: MIC value in mg/L for ciprofloxacin (numeric; n =
  5,360).

- `Cefixime`: MIC value in mg/L for cefixime (character, to preserve
  inequality prefixes such as `<0.016`).

- `Ceftriaxone`: MIC value in mg/L for ceftriaxone (numeric; n = 5,361).

## Source

ENA projects
[PRJEB9227](https://www.ebi.ac.uk/ena/browser/view/PRJEB9227),
[PRJEB34068](https://www.ebi.ac.uk/ena/browser/view/PRJEB34068),
[PRJEB58139](https://www.ebi.ac.uk/ena/browser/view/PRJEB58139).
