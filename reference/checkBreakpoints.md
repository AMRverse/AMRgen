# Check and Retrieve Breakpoints for an Antibiotic

This function checks the clinical breakpoints for a specified antibiotic
and species using the `getBreakpoints` function. It handles cases where
multiple breakpoint sites exist and uses the specified site or the one
with the highest susceptibility breakpoint.

## Usage

``` r
checkBreakpoints(
  species,
  guide = "EUCAST 2024",
  antibiotic,
  assay = "MIC",
  bp_site = NULL
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

  A character string indicating the antibiotic for which to check
  breakpoints (e.g., "Ciprofloxacin").

- assay:

  A character string specifying the assay type (either "MIC" or "Disk").
  Default is "MIC".

- bp_site:

  A character string specifying the breakpoint site to use (optional).
  If provided, the function uses this site; otherwise, if different
  breakpoints are specified for different sites it selects the one with
  the highest susceptibility breakpoint.

## Value

A list containing:

- `breakpoint_S`: The susceptibility breakpoint (e.g., MIC or disk
  size).

- `breakpoint_R`: The resistance breakpoint (e.g., MIC or disk size).

- `bp_standard`: The breakpoint site used, if multiple breakpoints are
  set in the guidelines.

## Details

The function first attempts to retrieve breakpoints using the
`getBreakpoints` function. If multiple breakpoint sites are found, it
handles the situation by:

- Using the specified site if it exists.

- Selecting the breakpoint with the highest susceptibility value if the
  specified site is not found.

- Returning a message about the selected site and breakpoint values.

## Examples

``` r
checkBreakpoints(species="Escherichia coli", guide="EUCAST 2024", 
                      antibiotic="Ciprofloxacin", assay="MIC")
#>   MIC breakpoints determined using AMR package: S <= 0.25 and R > 0.5
#>   NOTE: Multiple breakpoint entries, for different sites: Non-meningitis; Meningitis. Using the one with the highest S breakpoint (Non-meningitis).
#> $breakpoint_S
#> [1] 0.25
#> 
#> $breakpoint_R
#> [1] 0.5
#> 
#> $bp_standard
#> [1] "Non-meningitis"
#> 
```
