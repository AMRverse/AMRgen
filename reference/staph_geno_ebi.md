# S. aureus Example of Imported EBI Genotype Data

Phenotypes sourced from EBI AMR Portal using the
[download_ebi](https://AMRverse.github.io/AMRgen/reference/download_ebi.md)
function and imported to AMRgen phenotype table format.

## Usage

``` r
staph_geno_ebi
```

## Format

### `staph_geno_ebi` A data frame with 198344 rows and 34 columns representing genotyping results from EBI, imported into AMRgen format using [import_amrfp](https://AMRverse.github.io/AMRgen/reference/import_amrfp.md).

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

- ...: Additional data columns from EBI.

## Source

<https://www.ebi.ac.uk/amr>
