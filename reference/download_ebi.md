# Download antimicrobial genotype or phenotype data from the EBI AMR Portal

This function will retrieve genotype or phenotype data from the EBI AMR
Portal, via FTP. The portal uses AMRfinderplus to identify
AMR-associated genotypes, but the results are processed and not all
fields returned by AMRfinderplus are included. Optionally, the function
can also reformat the phenotype data for easy use with AMRgen functions
(using
[`import_ebi_ast_ftp()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast_ftp.md))
and re-interpret assay measures using the latest breakpoints/ECOFF.

## Usage

``` r
download_ebi(
  data = "phenotype",
  antibiotic = NULL,
  force_antibiotic = FALSE,
  genus = NULL,
  species = NULL,
  geno_subclass = NULL,
  geno_class = NULL,
  remove_dup = FALSE,
  release = NULL,
  reformat = FALSE,
  interpret_eucast = FALSE,
  interpret_clsi = FALSE,
  interpret_ecoff = FALSE
)
```

## Arguments

- data:

  String specifying the type of data to download, either "phenotype" or
  "genotype" (default "phenotype").

- antibiotic:

  (Optional) String (or vector of strings) specifying the antibiotic
  name/s to filter on (default NULL). Uses the AMR package to try to fix
  typos, and format to lower-case for EBI files. Not used if
  `data`="genotype" and `class` or \`subclass“ is specified.

- force_antibiotic:

  (Optional) Logical indicating whether to turn off parsing of
  antibiotic names and match exactly on the input strings (default
  FALSE).

- genus:

  (Optional) String specifying a bacterial genus to filter on (default
  NULL, will pull all taxa).

- species:

  (Optional) String specifying a bacterial species to filter on (default
  NULL, will pull all taxa). Not used if genus is specified.

- geno_subclass:

  (Optional) String specifying an antibiotic subclass to filter genotype
  data on (default NULL). Filter is based on string match, not identity,
  so e.g. subclass="TRIMETHOPRIM" will return all rows where the string
  "TRIMETHOPRIM" is included in the subclass field. Only used if
  `data`="genotype". Check NCBI Subclass list for valid terms
  (https://github.com/ncbi/amr/wiki/class-subclass).

- geno_class:

  (Optional) String specifying an antibiotic subclass to filter genotype
  data on (default NULL). Filter is based on string match, not identity,
  so e.g. class="TRIMETHOPRIM" will return all rows where the string
  "TRIMETHOPRIM" is included in the class field. Only used if
  `data`="genotype" and subclass is not specified. Check NCBI Class list
  for valid terms (https://github.com/ncbi/amr/wiki/class-subclass).

- remove_dup:

  (Optional) Logical specifying whether to clean up genotype data by
  removing duplicates for the same hit (default FALSE). Where a detected
  gene is associated with a subclass value that is actually a list of
  multiple drugs/classes, e.g. subclass="GENTAMICIN/TOBRAMYCIN", the EBI
  data table will have duplicate rows for the same gene hit, but with
  different values for the `antibiotic_name` and associated
  `antibiotic_ontology` and `antibiotic_ontology_link` annotation fields
  (e.g. one row each for gentamicin and tobramycin). To remove these
  duplicate rows (and the drug-specific annotation fields) and return
  only one row per hit (i.e. restoring AMRfinderplus output format), set
  this to TRUE.

- release:

  (Optional) String specifying the data release to download (default
  NULL, will pull latest release).

- reformat:

  (Optional) Logical specifying whether to reformat the downloaded
  phenotype data for easy use with downstream AMRgen function, using
  [`import_ebi_ast_ftp()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast_ftp.md).
  Default `FALSE`. This does things like format the antibiotic,
  measurement, and phenotype columns to AMR package classes. If set to
  `TRUE` you can also turn on re-interpreting MIC/disk data using latest
  EUCAST/CLSI breakpoints. No columns are removed from the downloaded
  data frame, but key fields are renamed, see documentation for
  [`format_ast()`](https://AMRverse.github.io/AMRgen/reference/format_ast.md).

- interpret_eucast:

  (Optional) Logical specifying whether to re-interpret the
  susceptibility phenotype (SIR) for each row based on the MIC or disk
  diffusion values, against EUCAST human breakpoints. These will be
  reported in a new column `pheno_eucast`, of class 'sir'. Only used
  when downloading phenotype data, with reformat set to `TRUE`.

- interpret_clsi:

  (Optional) Logical specifying whether to re-interpret the
  susceptibility phenotype (SIR) for each row based on the MIC or disk
  diffusion values, against CLSI human breakpoints. These will be
  reported in a new column `pheno_clsi`, of class 'sir'. Only used when
  downloading phenotype data, with reformat set to `TRUE`.

- interpret_ecoff:

  (Optional) Logical specifying whether to re-interpret the wildtype vs
  nonwildtype status for each row based on the MIC or disk diffusion
  values, against epidemiological cut-off (ECOFF) values. These will be
  reported in a new column `ecoff`, of class 'sir' and coded as 'R'
  (nonwildtype) or 'S' (wildtype). Only used when downloading phenotype
  data, with reformat set to `TRUE`.

## Value

A data frame containing EBI genotype data

## Details

See https://www.ebi.ac.uk/amr/about/ for more information on what is
available in the portal, and
https://github.com/ncbi/amr/wiki/class-subclass for valid class and
subclass terms.

Note the function downloads the full genotype or phenotype data table
before filtering on the provided parameters, so if you are having
trouble with drug/class names not matching then just run without
specifying any genus/species/antibiotic/class filters, to get the full
unfiltered table and explore the field values to filter manually to get
what you want.

## Examples
