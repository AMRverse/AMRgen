# Generate a Stacked Bar Plot of MIC Values Colored by Gene Symbol for Each Antibiotic

This function creates a stacked bar plot using `ggplot2`, where the
x-axis represents MIC (Minimum Inhibitory Concentration) or disk values,
the y-axis represents their frequency, and the bars are colored to
indicate whether the assay value is expressed as a range or not. It can
optionally be faceted on an additional categorical variable.

## Usage

``` r
assay_by_var(
  pheno_table,
  antibiotic,
  measure = "mic",
  var = "method",
  species = NULL,
  marker_free_strains = NULL,
  bp_site = NULL,
  cols = c(range = "maroon", value = "navy", `NA` = "grey")
)
```

## Arguments

- pheno_table:

  Phenotype table in standard format as per import_ast().

- antibiotic:

  Name of the antibiotic.

- measure:

  Field name containing the assay measurements to plot (default "mic").

- var:

  Field name containing a field to facet on (default "method").

- species:

  (optional) Name of species, so we can retrieve breakpoints to print at
  the top of the plot to help interpret it.

- marker_free_strains:

  (optional) Vector of sample names to select to get their own plot.
  Most useful for defining the set of strains with no known markers
  associated with the given antibiotic, so you can view the distribution
  of assay values for strains expected to be wildtype, which can help to
  identify issues with the assay.

- bp_site:

  (optional) Breakpoint site to retrieve (only relevant if also
  supplying species to retrieve breakpoints).

- cols:

  (optional) Manual colour scale to use for plot.

## Value

A list containing

- plot:

  Main plot with all samples that have assay data for the given
  antibiotic

- plot_nomarkers:

  Additional plot showing only those samples listed in
  `marker_free_strains`
