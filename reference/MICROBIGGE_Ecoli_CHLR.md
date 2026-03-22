# MicroBIGG-E for E.coli containing all chloramphenicol resistance genes

Data from NCBI Pathogen Detection Microbial Browser for Identification
of Genetic and Genomic Elements (MicroBIGG-E), containing all
chloramphenicol resistance gene hits in Escherichia coli genomes,
imported using
[`import_amrfp()`](https://amrgen.org/reference/import_amrfp.md).

## Usage

``` r
MICROBIGGE_Ecoli_CHLR
```

## Format

`MICROBIGGE_Ecoli_CHLR` A data frame with 95,776 rows and 27 columns:

- `id`: BioSample.

- `drug_agent`, `drug_class`: Antibiotic agent and class, determined by
  parsing AMRFinderPlus `subclass` field in the downloaded file.

- `gene`, `node`, `marker`: gene identifiers.

- `mutation`: mutation within gene, parsed into HGVS nomenclature format
  from `amr_element_symbol` field in the downloaded file.

- `% Coverage of reference`: % Coverage of reference.

- `% Identity to reference`: % Identity to reference.

- ...: Additional data columns from AMRFinderPlus \#' @source
  <https://www.ncbi.nlm.nih.gov/pathogens/microbigge/#chloramphenicol%20AND%20Escherichia>
