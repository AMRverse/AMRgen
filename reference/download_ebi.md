# Download antimicrobial genotype or phenotype data from the EBI AMR Portal

This function will retrieve genotype or phenotype data from the EBI AMR
Portal, via FTP. The portal uses AMRfinderplus to identify
AMR-associated genotypes, but the results are processed and not all
fields returned by AMRfinderplus are included. See
https://www.ebi.ac.uk/amr/about/#AMR-Genotypes for more information, and
https://github.com/ncbi/amr/wiki/class-subclass for valid class and
subclass terms. Optionally, the function can also reformat the phenotype
data for easy use with AMRgen functions (using
[`import_ebi_ast_ftp()`](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast_ftp.md))
and re-interpret assay measures using the latest breakpoints/ECOFF.

## Usage

``` r
download_ebi(
  data = "phenotype",
  antibiotic = NULL,
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

  (Optional) String specifying an antibiotic to to filter on (default
  NULL). Not used if class or subclass is specified.

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
  `data`="genotype". Check NCBI Subclass for valid terms.

- geno_class:

  (Optional) String specifying an antibiotic subclass to filter on
  (default NULL). Filter is based on string match, not identity, so e.g.
  class="TRIMETHOPRIM" will return all rows where the string
  "TRIMETHOPRIM" is included in the class field. Only used if
  `data`="genotype" and subclass is not specified. Check NCBI Class for
  valid terms.

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

## Examples

``` r
if (FALSE) { # \dontrun{
# download all phenotype data from EBI
pheno_ebi <- download_ebi()

# download phenotype data from Dec 2025 release, and filter to Salmonella
pheno_salmonella <- download_ebi(
    genus="Salmonella",
    user_release="2025-12"
)

# reformat downloaded phenotype data to simplify use with AMRgen functions
pheno_salmonella <- import_ebi_ast_ftp(pheno_salmonella)

# download phenotype data for Staphylococcus aureus and reformat 
# for use with AMRgen functions
pheno_staph <- download_ebi(
    species="Staphylococcus aureus",
    reformat=T
)

# download phenotype data for Klebsiella quasipneumoniae, reformat 
# for use with AMRgen functions, and re-interpret phenotypes
pheno_kquasi_reinterpreted <- download_ebi(
    species="Klebsiella quasipneumoniae",
    reformat=T,
    interpret_eucast = TRUE, 
    interpret_clsi = TRUE, 
    interpret_ecoff = TRUE
)

ebi_geno <- download_ebi(data="genotype")
    
geno_kpn_tmp <- download_ebi(
    data="genotype",
    species="Klebsiella pneumoniae",
    geno_subclass="TRIMETHOPRIM"
)
} # }
```
