# E. coli Ciprofloxacin MIC Distribution Example Data

Ciprofloxacin MIC distributions for E. coli, calculated from public data
and compared with the EUCAST reference distribution.

## Usage

``` r
ecoli_cip_vs_ref
```

## Format

An object of class `compare_eucast` with 32 rows and 3 columns. It
provides MIC distributions from EUCAST and public AST data extracted
from
[ecoli_ast](https://AMRverse.github.io/AMRgen/reference/ecoli_ast.md) in
the form of counts per value.

Columns include:

- `value`: MIC value.

- `user`: Count of samples with this MIC value, from the example data
  [ecoli_ast](https://AMRverse.github.io/AMRgen/reference/ecoli_ast.md).

- `eucast`: Count of samples with this MIC value, downloaded from EUCAST
  (Feb 2026).

## Source

<https://mic.eucast.org/>
