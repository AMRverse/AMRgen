# Example Resistance Gene Identifier (RGI) v6.0.6 Genotype Data from EuSCAPE project

Raw RGI v6.0.6 results file for Klebsiella pneumoniae from EuSCAPE
project, one AMR determinant per row. Includes only Perfect and Strict
hits Columns `Predicted_DNA`, `Predicted_Protein`, and
`CARD_Protein_Sequence` were removed to reduce file size.

## Usage

``` r
rgi_EuSCAPE_raw
```

## Format

`rgi_EuSCAPE_raw` A data frame with 59,403 rows and 25 columns:

- `ORF_ID`: Sample identifier

- ...: RGI results columns

## Source

ENA BioProject
[PRJEB10018](https://www.ebi.ac.uk/ena/browser/view/PRJEB10018). See
David *et al.* (2019) <https://doi.org/10.1038/s41564-019-0492-8>.
