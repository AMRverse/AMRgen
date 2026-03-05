# E. coli genotype data from Mills et al 2022

Genotyping data for isolates published in Mills et al, Genome Medicine
(2022) 14:147, generated using AMRfinderplus v3.12.8 and downloaded from
the
[AllTheBacteria](https://github.com/AllTheBacteria/AllTheBacteria/tree/main/reproducibility/All-samples/AMR/AMRFinderPlus)
project, and imported to AMRgen genotype table format. Corresponding MIC
data is available in `pheno_eco_2075`.

## Usage

``` r
geno_eco_2075
```

## Format

`geno_eco_2075` A data frame with 56064 rows and 24 columns representing
AMRfinderplus genotyping results for 2075 E. coli isolates.

Columns include:

- `id`: Sample identifier.

- `drug_agent`, `drug_class`: Antibiotic agent and class, determined by
  parsing AMRfinderplus `subclass` field in the downloaded file.

- `gene`, `node`, `marker`: gene symbol, parsed from
  `amr_element_symbol` field in the downloaded file.

- `mutation`: mutation within gene, parsed into HGVS nomenclature format
  from `amr_element_symbol` field in the downloaded file.

- `marker.label`: label for genotype marker, combining `gene` and
  `mutation` information (deletion variants represented as `"gene:-"`).

- ...: Additional data columns from AMRfinderplus

## Source

<https://github.com/AllTheBacteria/AllTheBacteria>
