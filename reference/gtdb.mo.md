# Get Microorganism from GTDB Species Name

Parse a character vector containing species names from GTDB output to
get a valid microorganism code.

## Usage

``` r
gtdb.mo(species)
```

## Arguments

- species:

  Name of species, coerced with
  [`AMR::as.mo()`](https://amr-for-r.org/reference/as.mo.html).

## Value

A character vector of valid microorganism codes as determined by
[`AMR::as.mo()`](https://amr-for-r.org/reference/as.mo.html).

## Examples

``` r
gtdb.mo("Escherichia_A coli_BC")
#> Class 'mo'
#> [1] B_ESCHR_COLI_COLI
```
