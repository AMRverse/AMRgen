# Import and process antimicrobial phenotype data from common sources

This function imports an antibiotic susceptibility testing datasets in
formats exported by EBI, NCBI, WHOnet and several automated AST
instruments (Vitek, Microscan, Sensititre, Phoenix). It optionally can
use the AMR package to interpret susceptibility phenotype (SIR) based on
EUCAST or CLSI guidelines (human breakpoints and/or ECOFF).

## Usage

``` r
import_pheno(
  input,
  format = "ebi",
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE,
  ...
)
```

## Arguments

- input:

  A string representing a dataframe, or a path to an input file,
  containing the phenotype data in a supported format. These files may
  be downloaded from public sources such as the [EBI AMR web
  browser](https://www.ebi.ac.uk/amr/data/?view=experiments), [EBI FTP
  site](ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/), or
  [NCBI browser](https://www.ncbi.nlm.nih.gov/pathogens/ast), or using
  the functions
  [`download_ebi()`](https://AMRverse.github.io/AMRgen/reference/download_ebi.md)
  or
  [`download_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/download_ncbi_ast.md);
  or the files may be exported from supported AST instruments.

- format:

  A string indicating the format of the data: `"ebi"` (default),
  `"ebi_web"`, `"ebi_ftp"`, `"ncbi"`, `"ncbi-biosample"`, `"vitek"`,
  `"microscan"`, `"phoenix"`, `"sensititre"`, or `"whonet"`. This
  determines which importer function the data is passed on to for
  processing (see below).

- interpret_eucast:

  A logical value (default is `FALSE`). If `TRUE`, the function will
  interpret the susceptibility phenotype (SIR) for each row based on the
  MIC or disk diffusion values, against EUCAST human breakpoints. These
  will be reported in a new column `pheno_eucast`, of class 'sir'.

- interpret_clsi:

  A logical value (default is `FALSE`). If `TRUE`, the function will
  interpret the susceptibility phenotype (SIR) for each row based on the
  MIC or disk diffusion values, against CLSI human breakpoints. These
  will be reported in a new column `pheno_clsi`, of class 'sir'.

- interpret_ecoff:

  A logical value (default is `FALSE`). If `TRUE`, the function will
  interpret the wildtype vs nonwildtype status for each row based on the
  MIC or disk diffusion values, against epidemiological cut-off (ECOFF)
  values. These will be reported in a new column `ecoff`, of class 'sir'
  and coded as `NWT` (nonwildtype) or `WT` (wildtype).

- ...:

  Format-specific arguments. See

  - `"ebi"` :
    [`import_ebi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast.md)

  - `"ebi_web"` :
    [`import_ebi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast.md)

  - `"ebi_ftp"`
    :[`import_ebi_ast_ftp()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast_ftp.md)

  - `"ncbi"` :
    [`import_ncbi_ast()`](https://AMRverse.github.io/AMRgen/reference/import_ncbi_ast.md)

  - `"ncbi_biosample"` :
    [`import_ncbi_biosample()`](https://AMRverse.github.io/AMRgen/reference/import_ncbi_biosample.md)

  - `"vitek"` :
    [`import_vitek_ast()`](https://AMRverse.github.io/AMRgen/reference/import_vitek_ast.md)

  - `"microscan"` :
    [`import_microscan_ast()`](https://AMRverse.github.io/AMRgen/reference/import_microscan_ast.md)

  - `"sensititre"` :
    [`import_sensititre_ast()`](https://AMRverse.github.io/AMRgen/reference/import_sensititre_ast.md)

  - `"phoenix"` :
    [`import_phoenix_ast()`](https://AMRverse.github.io/AMRgen/reference/import_phoenix_ast.md)

  - `"whonet"` :
    [`import_whonet_ast()`](https://AMRverse.github.io/AMRgen/reference/import_whonet_ast.md)

## Value

A data frame with the processed AST data, including additional columns:

- `id`: The biosample identifier (`character`).

- `spp_pheno`: The species phenotype, formatted using the
  [`AMR::as.mo()`](https://amr-for-r.org/reference/as.mo.html) function
  (class `mo`).

- `drug_agent`: The antibiotic used in the test, formatted using the
  [`AMR::as.ab()`](https://amr-for-r.org/reference/as.ab.html) function
  (class `ab`).

- `mic`: The minimum inhibitory concentration (MIC) value, formatted
  using the
  [`AMR::as.mic()`](https://amr-for-r.org/reference/as.mic.html)
  function (class `mic`).

- `disk`: The disk diffusion measurement (in mm), formatted using the
  [`AMR::as.disk()`](https://amr-for-r.org/reference/as.disk.html)
  function (class `disk`).

- `method`: The AST method (e.g., `"broth dilution"`,
  `"disk diffusion"`, `"Etest"`, `"agar dilution"`). Expected values are
  based on the NCBI/EBI antibiogram specification (`character`).

- `platform`: The AST platform/instrument (e.g., `"Vitek"`, `"Phoenix"`,
  `"Sensititre"`) (`character`).

- `guideline`: The AST standard recorded in the input file as being used
  for the AST assay (`character`).

- `pheno_eucast`: The phenotype newly interpreted against EUCAST human
  breakpoint standards (as `S/I/R`), based on the MIC or disk diffusion
  data (class `sir`).

- `pheno_clsi`: The phenotype newly interpreted against CLSI human
  breakpoint standards (as `S/I/R`), based on the MIC or disk diffusion
  data (class `sir`).

- `ecoff`: The phenotype newly interpreted against the ECOFF (as
  `WT/NWT`), based on the MIC or disk diffusion data (class `sir`).

- `pheno_provided`: The original phenotype interpretation provided in
  the input file, formatted using
  [`AMR::as.sir()`](https://amr-for-r.org/reference/as.sir.html) (class
  `sir`).

- `source`: The source of each data point (from the publications or
  bioproject field in the input file, or replaced with a single value
  passed in as the `source` parameter) (`character`).

## Examples

``` r
if (FALSE) { # \dontrun{
# import NCBI data without re-interpreting resistance
pheno <- import_pheno(staph_ast_ncbi_cloud_raw, format = "ncbi")

# import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
pheno <- import_ast(ecoli_ast_raw, format = "ncbi", interpret_eucast = TRUE, interpret_ecoff = TRUE)
head(pheno)

# download EBI phenotype data for Klebsiella quasipneumoniae
kquasi_raw_ebi <- download_ebi(
  species = "Klebsiella quasipneumoniae"
)
# import the data and interpret against ecoff
pheno <- import_pheno(kquasi_raw_ebi,
  format = "ebi_ftp",
  interpret_ecoff = TRUE
)
# Download NCBI phenotype data for Klebsiella quasipneumoniae
kquasi_raw_ncbi <- download_ncbi_ast(
  "Klebsiella quasipneumoniae"
)
# import the data and interpret against EUCAST breakpoints
ast <- import_pheno(kquasi_raw_ncbi, 
  format = "ncbi_biosample",
  interpret_eucast = T)

# import Vitek data from file, with default parameters
pheno <- import_pheno("vitek_export.tsv",
  format = "vitek"
)

# import Vitek data from file
# specify guideline that was used, remove dates, ignore expertized calls
pheno <- import_pheno("vitek_export.tsv",
  format = "vitek",
  instrument_guideline = "EUCAST 2025",
  use_expertized = FALSE,
  include_dates = FALSE
)
} # }
```
