# Import GTDB Output

Import GTDB output (from file or data frame) and parse the species name
into a microorganism recognised by AMR package. AMR mo code and species
name will be added to the data frame.

## Usage

``` r
import_gtdb(file = NULL, tbl = NULL, species_column = "Species")
```

## Arguments

- file:

  File path to GTDB output file (tab-separated). If not given, `tbl`
  must be provided.

- tbl:

  Data frame containing GTDB output. If not given, `file` must be
  provided.

- species_column:

  Name of column containing the species call (default "Species").

## Value

A data frame containing GTDB output with AMR-parsed microorganism code
(`gtdb.mo`) and species name (`gtdb.species`) appended.

## Examples

``` r
if (FALSE) { # \dontrun{
import_gtdb(tbl = data.frame(Species = c(
  "Pseudomonas_E piscis",
  "Haemophilus_D parainfluenzae_A",
  "Acinetobacter calcoaceticus_C"
)))
} # }
```
