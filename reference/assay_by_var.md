# Generate a Stacked Bar Plot of Assay Values Colored by a Variable

This function creates a stacked bar plot using `ggplot2`, where the
x-axis represents MIC (Minimum Inhibitory Concentration) or disk values,
the y-axis indicates their frequency, and the bars are colored by a
variable (by default, colours indicate whether the assay value is
expressed as a range or not). Plots can optionally be faceted on an
additional categorical variable. If breakpoints are provided, or species
and drug are provided so we can extract EUCAST breakpoints, vertical
lines indicating the S/R breakpoints and ECOFF will be added to the
plot.

## Usage

``` r
assay_by_var(
  pheno_table,
  antibiotic,
  measure = "mic",
  facet_var = NULL,
  species = NULL,
  bp_site = NULL,
  bp_S = NULL,
  bp_R = NULL,
  ecoff = NULL,
  guideline = "EUCAST 2025",
  marker_free_strains = NULL,
  colour_by = NULL,
  cols = NULL
)
```

## Arguments

- pheno_table:

  Phenotype table in standard format as per import_ast().

- antibiotic:

  Name of the antibiotic.

- measure:

  Field name containing the assay measurements to plot (default "mic").

- facet_var:

  (optional) Field name containing a field to facet on (default NULL).

- species:

  (optional) Name of species, so we can retrieve breakpoints to print at
  the top of the plot to help interpret it.

- bp_site:

  (optional) Breakpoint site to retrieve (only relevant if also
  supplying `species` to retrieve breakpoints, and not supplying
  breakpoints via `bp_S`, `bp_R`, `ecoff`).

- bp_S:

  (optional) S breakpoint

- bp_R:

  (optional) R breakpoint

- ecoff:

  (optional) ECOFF breakpoint

- guideline:

  (optional) Guideline to use when looking up breakpoints (default
  'EUCAST 2025')

- marker_free_strains:

  (optional) Vector of sample names to select to get their own plot.
  Most useful for defining the set of strains with no known markers
  associated with the given antibiotic, so you can view the distribution
  of assay values for strains expected to be wildtype, which can help to
  identify issues with the assay.

- colour_by:

  (optional) Field name containing a field to colour bars by (default
  NULL, which will colour each bar to indicate whether the value is
  expressed as a range or not)

- cols:

  (optional) Manual colour scale to use for plot. If NULL, `colour_by`
  variable is of class 'sir', bars will by default be coloured using
  standard SIR colours.

## Value

A list containing

- plot:

  Main plot with all samples that have assay data for the given
  antibiotic

- plot_nomarkers:

  Additional plot showing only those samples listed in
  `marker_free_strains`

## Examples

``` r
# plot MIC distribution, highlighting values expressed as ranges
assay_by_var(pheno_table=ecoli_ast, antibiotic="Ciprofloxacin", 
                measure="mic")
#> $plot_nomarkers
#> NULL
#> 
#> $plot

#> 

# colour by SIR interpretion recorded in column 'pheno_clsi'
assay_by_var(pheno_table=ecoli_ast, antibiotic="Ciprofloxacin", 
                measure="mic", colour_by = "pheno_clsi")
#> $plot_nomarkers
#> NULL
#> 
#> $plot

#> 

# look up ECOFF and CLSI breakpoints and annotate these on the plot
assay_by_var(pheno_table=ecoli_ast, antibiotic="Ciprofloxacin", 
                measure="mic", colour_by = "pheno_clsi", 
                species="E. coli", guideline="CLSI 2025")
#> Error in executing command: object of type 'builtin' is not subsettable
#>   MIC breakpoints determined using AMR package: S <= 0.25 and R > 1
#> $plot_nomarkers
#> NULL
#> 
#> $plot
#> Warning: Removed 26 rows containing missing values or values outside the scale range
#> (`geom_vline()`).

#> 

# facet by method
assay_by_var(pheno_table=ecoli_ast, antibiotic="Ciprofloxacin", 
                measure="mic", colour_by = "pheno_clsi", 
                species="E. coli", guideline="CLSI 2025", 
                facet_var ="method")
#> Error in executing command: object of type 'builtin' is not subsettable
#>   MIC breakpoints determined using AMR package: S <= 0.25 and R > 1
#> $plot_nomarkers
#> NULL
#> 
#> $plot
#> Warning: Removed 208 rows containing missing values or values outside the scale range
#> (`geom_vline()`).

#> 
```
