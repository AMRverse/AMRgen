# Import/Export BioSample Antibiograms

Output phenotype data to [NCBI BioSample
Antibiograms](https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/).

## Usage

``` r
export_ncbi_biosample(data, file, overwrite = FALSE)
```

## Arguments

- data:

  Data set containing SIR results.

- file:

  File path to which the data will be exported.

- overwrite:

  Logical indicating whether to overwrite an existing file.

## Value

A tab-delimited UTF-8 text file is written to disk in the format
required by NCBI BioSample Antibiograms.
