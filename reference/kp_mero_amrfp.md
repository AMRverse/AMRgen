# Example AMRFinderPlus Genotype Data from EuSCAPE project

AMRFinderPlus results file for Klebsiella pneumoniae from EuSCAPE
project, one AMR determinant per row, downloaded from the EBI AMR portal
using [`download_ebi()`](https://amrgen.org/reference/download_ebi.md)
and imported using
[`import_geno()`](https://amrgen.org/reference/import_geno.md).

## Usage

``` r
kp_mero_amrfp
```

## Format

`kp_mero_amrfp` A data frame with 32,385 rows and 34 columns:

- `id`: BioSample.

- `drug`, `drug_class`: Antibiotic agent and class, determined by
  parsing AMRFinderPlus `subclass` field in the downloaded file.

- `gene`, `node`, `marker`: gene identifiers.

- `mutation`: mutation within gene, parsed into HGVS nomenclature format
  from `amr_element_symbol` field in the downloaded file.

- `% Coverage of reference`: % Coverage of reference.

- `% Identity to reference`: % Identity to reference.

- ...: Additional data columns from AMRFinderPlus

## Source

[EBI AMR Portal](https://www.ebi.ac.uk/amr). See David *et al.* (2019)
<https://doi.org/10.1038/s41564-019-0492-8>.
