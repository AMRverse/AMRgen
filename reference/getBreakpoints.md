# Get Clinical Breakpoints for an Antibiotic

This function retrieves the clinical breakpoints for a given species,
antibiotic, and guideline, from the AMR package. It attempts to find the
breakpoints at various taxonomic levels (species, genus, family, order)
if no direct match is found.

## Usage

``` r
getBreakpoints(
  species,
  guide = "EUCAST 2024",
  antibiotic,
  type_filter = "human"
)
```

## Arguments

- species:

  A character string representing the species of interest (e.g.,
  "Escherichia coli").

- guide:

  A character string indicating the guideline for breakpoints (default,
  "EUCAST 2024").

- antibiotic:

  A character string indicating the antibiotic for which to retrieve
  breakpoints (e.g., "Ciprofloxacin").

- type_filter:

  A character string indicating the type of breakpoints to retrieve
  (e.g., "human"). Default is "human" which returns human clinical
  breakpoints, change to "ECOFF" to get the epidemiological cutoff.

## Value

A data frame containing the clinical breakpoints for the specified
species, antibiotic, and guideline. If no exact match is found, the
function attempts to retrieve breakpoints at the genus, family, and
order levels.

## Examples

``` r
getBreakpoints("Escherichia coli", "EUCAST 2024", "Ciprofloxacin")
#> # A tibble: 3 × 14
#>   guideline   type  host  method site  mo               rank_index ab   ref_tbl 
#>   <chr>       <chr> <chr> <chr>  <chr> <mo>                  <dbl> <ab> <chr>   
#> 1 EUCAST 2024 human human DISK   Non-… B_[ORD]_ENTRBCTR          5 CIP  Enterob…
#> 2 EUCAST 2024 human human MIC    Non-… B_[ORD]_ENTRBCTR          5 CIP  Enterob…
#> 3 EUCAST 2024 human human MIC    Meni… B_[ORD]_ENTRBCTR          5 CIP  Enterob…
#> # ℹ 5 more variables: disk_dose <chr>, breakpoint_S <dbl>, breakpoint_R <dbl>,
#> #   uti <lgl>, is_SDD <lgl>
```
