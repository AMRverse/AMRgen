# Meropenem Phenotype Data from EuSCAPE project

Meropenem phenotype data for Klebsiella pneumoniae from EuSCAPE project,
one sample per row, downloaded from the EBI AMR portal using
[`download_ebi()`](https://amrgen.org/reference/download_ebi.md) and
imported using
[`import_pheno()`](https://amrgen.org/reference/import_pheno.md).

## Usage

``` r
kp_mero_euscape
```

## Format

`kp_mero_euscape` A data frame with 1,490 rows and 43 columns:

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

- ...: Additional data columns from EBI AMR Portal

## Source

[EBI AMR Portal](https://www.ebi.ac.uk/amr). See David *et al.* (2019)
<https://doi.org/10.1038/s41564-019-0492-8>.
