# Import phenotype predictions from AMRrules output

This imports antimicrobial phenotype predictions (S/I/R and WT/NWT)
generated from genotype data using [AMRrules](https://www.amrrules.org).

## Usage

``` r
import_amrrules_predictions(
  input,
  sample_col = "sample",
  species_col = "organism",
  ab_col = "drug",
  sir_col = "clinical category",
  ecoff_col = "phenotype",
  method = "genotyping",
  platform = "AMRfinderplus + AMRrules"
)
```

## Arguments

- input:

  A string representing a dataframe, or a path to an input file,
  containing the AMRrules output "genome_summary" file, which should be
  a long-form TSV file with one row per sample and drug.

- sample_col:

  (optional, default `"sample"`) String indicating the name of the input
  data column that provides the sample name.

- species_col:

  (optional, default `"organism"`) String indicating the name of the
  input data column that provides a species name. If provided, this
  column will be converted to micro-organism class `mo` via
  [`AMR::as.mo()`](https://amr-for-r.org/reference/as.mo.html). If the
  `rename_cols` parameter is set to `TRUE`, this column will also be
  renamed as `spp_pheno`. If interpretation is switched on, this column
  will be used to identify the appropriate breakpoints for
  interpretation of each row in the data table.

- ab_col:

  (optional, default `"drug"`) String indicating the name of the input
  data column that provides a drug name.

- sir_col:

  (optional, default `"clinical category"`) String indicating the name
  of the input data column that indicates the S/I/R prediction.

- ecoff_col:

  (optional, default `"phenotype"`) String indicating the name of the
  input data column that indicates the WT/NWT prediction.

- method:

  (optional, default `"genotyping"`) String indicating the value to
  record in a new `method` field added to the output table.

- platform:

  (optional, default `"AMRfinderplus + AMRrules"`) String indicating the
  value to record in a new `platform` field added to the output table.

## Value

A data frame with the processed AST data, including additional columns:

## Examples

``` r
if (FALSE) { # \dontrun{
# import and process AST data from EBI, write formatted data to file for later use
predictions <- import_amrrules_predictions("Ecoli_genome_summary.tsv")
} # }
```
