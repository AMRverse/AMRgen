# Download antimicrobial genotype or phenotype data from the EBI AMR Portal

This function will retrieve genotype or phenotype data from the EBI AMR
Portal FTP site. The portal uses AMRfinderplus to identify
AMR-associated genotypes, but the results are processed and not all
fields returned by AMRfinderplus are included. Optionally, the function
can also reformat the phenotype data for easy use with AMRgen functions
(using
[import_ebi_ast_ftp](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast_ftp.md))
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
  `data`="genotype" and `class` or `subclass` is specified.

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
  `data`="genotype". Check [NCBI AMR Class-Subclass
  Reference](https://github.com/ncbi/amr/wiki/class-subclass) for valid
  terms.

- geno_class:

  (Optional) String specifying an antibiotic subclass to filter genotype
  data on (default NULL). Filter is based on string match, not identity,
  so e.g. class="TRIMETHOPRIM" will return all rows where the string
  "TRIMETHOPRIM" is included in the class field. Only used if
  `data`="genotype" and subclass is not specified. Check [NCBI AMR
  Class-Subclass
  Reference](https://github.com/ncbi/amr/wiki/class-subclass) for valid
  terms.

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
  this to `TRUE`.

- release:

  (Optional) String specifying the data release to download (default
  NULL, will pull latest release).

- reformat:

  (Optional) Logical specifying whether to reformat the downloaded data
  for easy use with downstream AMRgen functions, using
  [import_ebi_ast_ftp](https://AMRverse.github.io/AMRgen/reference/import_ebi_ast_ftp.md)
  (phenotypes) or
  [import_amrfp_ebi_ftp](https://AMRverse.github.io/AMRgen/reference/import_amrfp_ebi_ftp.md).
  Default `FALSE`. This does things like format the antibiotic,
  measurement, and phenotype columns to AMR package classes. If set to
  `TRUE` you can also turn on re-interpreting MIC/disk data using latest
  EUCAST/CLSI breakpoints (when `data`="phenotype"). No columns are
  removed from the downloaded data frame, but key fields are renamed,
  see documentation for
  [format_ast](https://AMRverse.github.io/AMRgen/reference/format_ast.md)
  and
  [import_amrfp_ebi_ftp](https://AMRverse.github.io/AMRgen/reference/import_amrfp_ebi_ftp.md).

- interpret_eucast:

  (Optional) Logical specifying whether to re-interpret the
  susceptibility phenotype (SIR) for each row based on the MIC or disk
  diffusion values, against EUCAST human breakpoints. These will be
  reported in a new column `pheno_eucast`, of class `sir`. Only used
  when downloading phenotype data, with reformat set to `TRUE`.

- interpret_clsi:

  (Optional) Logical specifying whether to re-interpret the
  susceptibility phenotype (SIR) for each row based on the MIC or disk
  diffusion values, against CLSI human breakpoints. These will be
  reported in a new column `pheno_clsi`, of class `sir`. Only used when
  downloading phenotype data, with reformat set to `TRUE`.

- interpret_ecoff:

  (Optional) Logical specifying whether to re-interpret the wildtype vs
  nonwildtype status for each row based on the MIC or disk diffusion
  values, against epidemiological cut-off (ECOFF) values. These will be
  reported in a new column `ecoff`, of class `sir` and coded as `R`
  (nonwildtype) or `S` (wildtype). Only used when downloading phenotype
  data, with reformat set to `TRUE`.

## Value

A data frame containing the phenotype or genotype data retrieved from
EBI, optionally reformatted to AMRgen standard formats and classes.

## Details

See <https://www.ebi.ac.uk/amr/about/> for more information on what is
available in the portal, and
<https://github.com/ncbi/amr/wiki/class-subclass> for valid class and
subclass terms.

Note the function downloads the full genotype or phenotype data table
before filtering on the provided parameters, so if you are having
trouble with drug/class names not matching then just run without
specifying any genus/species/antibiotic/class filters, to get the full
unfiltered table and explore the field values to filter manually to get
what you want.

## Examples

``` r
if (FALSE) { # \dontrun{
# download all phenotype data from EBI
pheno_ebi <- download_ebi()

# download phenotype data from Dec 2025 release, and filter to Salmonella
pheno_salmonella <- download_ebi(
  genus = "Salmonella",
  release = "2025-12"
)

# reformat downloaded phenotype data to simplify use with AMRgen functions
pheno_salmonella <- import_ebi_ast_ftp(pheno_salmonella)


# download phenotype data for Salmonella, filter to ampicillin and ciprofloxacin
pheno_salmonella <- download_ebi(
  genus = "Salmonella",
  antibiotic = c("ampicillin", "Cipro")
)

# download phenotype data for Staphylococcus aureus and reformat
# for use with AMRgen functions
pheno_staph <- download_ebi(
  species = "Staphylococcus aureus",
  reformat = T
)

# download phenotype data for Klebsiella quasipneumoniae, reformat
# for use with AMRgen functions, and re-interpret phenotypes
pheno_kquasi_reinterpreted <- download_ebi(
  species = "Klebsiella quasipneumoniae",
  reformat = T,
  interpret_eucast = TRUE,
  interpret_clsi = TRUE,
  interpret_ecoff = TRUE
)

# download all available genotype data
ebi_geno <- download_ebi(data = "genotype")

# download genotype data for Klebsiella pneumoniae, and filter to 
# markers assigned to NCBI class 'TRIMETHOPRIM'
geno_kpn_tmp <- download_ebi(
  data = "genotype",
  species = "Klebsiella pneumoniae",
  geno_subclass = "TRIMETHOPRIM"
)

# download genotype data for Klebsiella pneumoniae, and filter to 
# markers assigned to NCBI class 'TRIMETHOPRIM' or 'QUINOLONE'
geno_kpn_tmp <- download_ebi(
  data = "genotype",
  species = "Klebsiella pneumoniae",
  geno_subclass = c("TRIMETHOPRIM", "QUINOLONE")
)

# download genotype data for Klebsiella pneumoniae, and filter to 
# markers associated with CARD drug term 'trimethoprim'
geno_kpn_tmp <- download_ebi(
  data = "genotype",
  species = "Klebsiella pneumoniae",
  antibiotic = "trimethoprim"
)

} # }
```
