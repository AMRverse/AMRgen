# NCBI Subclass mapping to drug class

Mapping of NCBI refgene / AMRfinderplus Subclass terms that are not
present in the AMR package as drug class terms. Used internally when
importing AMRfinderplus results into AMRgen genotype table format.

## Usage

``` r
amrfp_drugs_table
```

## Format

`amrfp_drugs_table` A data frame with 21 rows and 2 columns.

Columns include:

- `AMRFP_Subclass`: NCBI term

- `drug_class`: Name of drug class as it should appear in imported
  genotype table.
