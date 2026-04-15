# Clindamycin MIC data for 5914 *Staphylococcus aureus* isolates

Used to investigate genetic determinants of clindamycin resistance.
Downloaded from NCBI and EBI.

## Usage

``` r
pheno_CLI_public
```

## Format

`pheno_CLI_public` A data frame with 5914 rows and 34 columns:

- `id`: Sample identifier, imported from the `#BioSample` column in the
  raw input.

- `drug`: Antibiotic code, interpreted from `Antibiotic` using `as.ab`,
  used to interpret `ecoff` and `pheno` columns.

- `mic`: Minimum inhibitory concentration, formatted using `as.mic`,
  used to interpret `ecoff` and `pheno` columns.

- `disk`: Disk diffusion zone, formatted using `as.disk`, used to
  interpret `ecoff` and `pheno` columns.

- `pheno_eucast`: S/I/R classification according to EUCAST, interpreted
  using `as.sir`.

- `ecoff`: WT/NWT classification, interpreted using `as.sir`.

- `guideline`: Interpretation guidelines used to interpret `ecoff` and
  `pheno` columns.

- `method`: Test method, one of: "E-test', "Unknown', "agar dilution',
  "broth dilution', "disk diffusion'.

- `platform`: Testing platform, one of "BD Phoenix', "Microscan',
  "Phoenix', "Sensititre', "Vitek'.

- `pheno_provided`: S/I/R interpretation as provided in the raw input.

- `spp_pheno`: Species identifier, interpreted from `Scientific name`
  using `as.mo`, used to interpret `ecoff` and `pheno` columns.

## Source

<https://www.ncbi.nlm.nih.gov/pathogens/ast>

<https://www.ebi.ac.uk/amr>
