# Convert single-letter amino acid code(s) to three-letter code(s)

This function takes a single-letter amino acid code, a vector of
single-letter codes, or a string representing a sequence of
single-letter codes. It returns the corresponding three-letter code(s),
concatenated directly for sequences.

## Usage

``` r
convert_aa_code(input_code)
```

## Arguments

- input_code:

  A character string (e.g., "A", "MAG") or a vector of character strings
  (e.g., c("A", "G", "C")).

## Value

A character string (or vector of strings) with the three-letter amino
acid code(s). For multi-character input strings, a single concatenated
string is returned. Returns NA for individual unmatched codes.

## Examples

``` r
# Single character input
convert_aa_code("A")
#> [1] "Ala"

# Vector of single characters
convert_aa_code(c("M", "A", "G", "Z")) # Z will be NA
#> [1] "Met" "Ala" "Gly" NA   

# Multi-character sequence input
convert_aa_code("MAG")
#> [1] "MetAlaGly"
convert_aa_code("MAGL")
#> [1] "MetAlaGlyLeu"
```
