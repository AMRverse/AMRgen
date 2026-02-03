# Convert mutation string based on method

This function takes a mutation string (e.g., "gene_REF123ALT") and a
mutation method, then extracts and converts parts of the mutation
string. Specifically designed for use within
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).

## Usage

``` r
convert_mutation(symbol_col, method_col)
```

## Arguments

- symbol_col:

  A character vector representing the 'Gene symbol' column (the mutation
  string).

- method_col:

  A character vector representing the 'Method' column.

## Value

A character vector containing the formatted mutation strings (e.g.,
"Ala123Trp") or NA if not applicable/match.
