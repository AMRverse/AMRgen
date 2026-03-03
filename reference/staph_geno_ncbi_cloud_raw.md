# S. aureus Example of Raw Genotype Data Downloaded from NCBI via Google Cloud BigQuery

AMRfinderplus genotypes sourced from NCBI via
[query_ncbi_bq_ast](https://AMRverse.github.io/AMRgen/reference/query_ncbi_bq_ast.md)
function, without reformating.

## Usage

``` r
staph_geno_ncbi_cloud_raw
```

## Format

`staph_geno_ncbi_cloud_raw` A data frame with 4064 rows and 9 columns
representing all Staphylococcus aureus genotyping results for markers
associated with class aminoglycoside or tetracycline.

Columns include:

- `biosample_acc`: Sample identifier.

- `scientific_name`: Organism name.

- `Gene symbol`, `Class`, `Subclass`, `Element type`, `Element subtype`,
  `Method`, `Hierarchy_node`: Key results fields from AMRfinderplus.

## Source

<https://www.ncbi.nlm.nih.gov/pathogens/microbigge/>
