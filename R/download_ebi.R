# ===================================================================== #
#  Licensed as GPL-v3.0.                                                #
#                                                                       #
#  Developed as part of the AMRverse (https://github.com/AMRverse):     #
#  https://github.com/AMRverse/AMRgen                                   #
#                                                                       #
#  We created this package for both routine data analysis and academic  #
#  research and it was publicly released in the hope that it will be    #
#  useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                       #
#  This R package is free software; you can freely use and distribute   #
#  it for both personal and commercial purposes under the terms of the  #
#  GNU General Public License version 3.0 (GNU GPL-3), as published by  #
#  the Free Software Foundation.                                        #
# ===================================================================== #
#' Download antimicrobial genotype or phenotype data from the EBI AMR Portal
#'
#' This function will retrieve genotype or phenotype data from the EBI AMR Portal FTP site. The portal uses AMRfinderplus to identify AMR-associated genotypes, but the results are processed and not all fields returned by AMRfinderplus are included.
#' Optionally, the function can also reformat the phenotype data for easy use with AMRgen functions (using [import_ebi_ast_ftp]) and re-interpret assay measures using the latest breakpoints/ECOFF.
#'
#' See <https://www.ebi.ac.uk/amr/about/> for more information on what is available in the portal, and <https://github.com/ncbi/amr/wiki/class-subclass> for valid class and subclass terms.
#'
#' Note the function downloads the full genotype or phenotype data table before filtering on the provided parameters, so if you are having trouble with drug/class names not matching then just run without specifying any genus/species/antibiotic/class filters, to get the full unfiltered table and explore the field values to filter manually to get what you want.
#'
#' @param data String specifying the type of data to download, either "phenotype" or "genotype" (default "phenotype").
#' @param antibiotic (Optional) String (or vector of strings) specifying the antibiotic name/s to filter on (default NULL). Uses the AMR package to try to fix typos, and format to lower-case for EBI files. Not used if `data`="genotype" and `class` or `subclass` is specified.
#' @param force_antibiotic (Optional) Logical indicating whether to turn off parsing of antibiotic names and match exactly on the input strings (default `FALSE`).
#' @param genus (Optional) String specifying a bacterial genus to filter on (default NULL, will pull all taxa).
#' @param species (Optional) String specifying a bacterial species to filter on (default NULL, will pull all taxa). Not used if genus is specified.
#' @param geno_subclass (Optional) String specifying an antibiotic subclass to filter genotype data on (default NULL). Filter is based on string match, not identity, so e.g. subclass="TRIMETHOPRIM" will return all rows where the string "TRIMETHOPRIM" is included in the subclass field. Only used if `data`="genotype". Check [NCBI AMR Class-Subclass Reference](https://github.com/ncbi/amr/wiki/class-subclass) for valid terms.
#' @param geno_class (Optional) String specifying an antibiotic subclass to filter genotype data on (default NULL). Filter is based on string match, not identity, so e.g. class="TRIMETHOPRIM" will return all rows where the string "TRIMETHOPRIM" is included in the class field. Only used if `data`="genotype" and subclass is not specified. Check [NCBI AMR Class-Subclass Reference](https://github.com/ncbi/amr/wiki/class-subclass) for valid terms.
#' @param remove_dup (Optional) Logical specifying whether to clean up genotype data by removing duplicates for the same hit (default FALSE). Where a detected gene is associated with a subclass value that is actually a list of multiple drugs/classes, e.g. subclass="GENTAMICIN/TOBRAMYCIN", the EBI data table will have duplicate rows for the same gene hit, but with different values for the `antibiotic_name` and associated `antibiotic_ontology` and `antibiotic_ontology_link` annotation fields (e.g. one row each for gentamicin and tobramycin). To remove these duplicate rows (and the drug-specific annotation fields) and return only one row per hit (i.e. restoring AMRfinderplus output format), set this to `TRUE`.
#' @param reformat (Optional) Logical specifying whether to reformat the downloaded data for easy use with downstream AMRgen functions, using [import_ebi_ast_ftp] (phenotypes) or [import_amrfp_ebi_ftp]. Default `FALSE`. This does things like format the antibiotic, measurement, and phenotype columns to AMR package classes. If set to `TRUE` you can also turn on re-interpreting MIC/disk data using latest EUCAST/CLSI breakpoints (when `data`="phenotype"). No columns are removed from the downloaded data frame, but key fields are renamed, see documentation for [format_ast] and [import_amrfp_ebi_ftp].
#' @param interpret_eucast (Optional) Logical specifying whether to re-interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`. Only used when downloading phenotype data, with reformat set to `TRUE`.
#' @param interpret_clsi (Optional) Logical specifying whether to re-interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`. Only used when downloading phenotype data, with reformat set to `TRUE`.
#' @param interpret_ecoff (Optional) Logical specifying whether to re-interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `R` (nonwildtype) or `S` (wildtype). Only used when downloading phenotype data, with reformat set to `TRUE`.
#' @param release (Optional) String specifying the data release to download (default NULL, will pull latest release).
#' @importFrom arrow read_parquet
#' @importFrom RCurl getBinaryURL getURL
#' @return A data frame containing the phenotype or genotype data retrieved from EBI, optionally reformatted to AMRgen standard formats and classes.
#' @export
#' @examples
#' \dontrun{
#' # download all phenotype data from EBI
#' pheno_ebi <- download_ebi()
#'
#' # download phenotype data from Dec 2025 release, and filter to Salmonella
#' pheno_salmonella <- download_ebi(
#'   genus = "Salmonella",
#'   release = "2025-12"
#' )
#'
#' # reformat downloaded phenotype data to simplify use with AMRgen functions
#' pheno_salmonella <- import_ebi_ast_ftp(pheno_salmonella)
#'
#'
#' # download phenotype data for Salmonella, filter to ampicillin and ciprofloxacin
#' pheno_salmonella <- download_ebi(
#'   genus = "Salmonella",
#'   antibiotic = c("ampicillin", "Cipro")
#' )
#'
#' # download phenotype data for Staphylococcus aureus and reformat
#' # for use with AMRgen functions
#' pheno_staph <- download_ebi(
#'   species = "Staphylococcus aureus",
#'   reformat = T
#' )
#'
#' # download phenotype data for Klebsiella quasipneumoniae, reformat
#' # for use with AMRgen functions, and re-interpret phenotypes
#' pheno_kquasi_reinterpreted <- download_ebi(
#'   species = "Klebsiella quasipneumoniae",
#'   reformat = T,
#'   interpret_eucast = TRUE,
#'   interpret_clsi = TRUE,
#'   interpret_ecoff = TRUE
#' )
#'
#' # download all available genotype data
#' ebi_geno <- download_ebi(data = "genotype")
#'
#' # download genotype data for Klebsiella pneumoniae, and filter to
#' # markers assigned to NCBI class 'TRIMETHOPRIM'
#' geno_kpn_tmp <- download_ebi(
#'   data = "genotype",
#'   species = "Klebsiella pneumoniae",
#'   geno_subclass = "TRIMETHOPRIM"
#' )
#'
#' # download genotype data for Klebsiella pneumoniae, and filter to
#' # markers assigned to NCBI class 'TRIMETHOPRIM' or 'QUINOLONE'
#' geno_kpn_tmp <- download_ebi(
#'   data = "genotype",
#'   species = "Klebsiella pneumoniae",
#'   geno_subclass = c("TRIMETHOPRIM", "QUINOLONE")
#' )
#'
#' # download genotype data for Klebsiella pneumoniae, and filter to
#' # markers associated with CARD drug term 'trimethoprim'
#' geno_kpn_tmp <- download_ebi(
#'   data = "genotype",
#'   species = "Klebsiella pneumoniae",
#'   antibiotic = "trimethoprim"
#' )
#' }
download_ebi <- function(data = "phenotype",
                         antibiotic = NULL, force_antibiotic = FALSE,
                         genus = NULL, species = NULL,
                         geno_subclass = NULL, geno_class = NULL,
                         remove_dup = FALSE,
                         release = NULL, reformat = FALSE,
                         interpret_eucast = FALSE,
                         interpret_clsi = FALSE,
                         interpret_ecoff = FALSE) {
  # EBI source url
  ebi_url <- "ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/"

  filename <- paste0(data, ".parquet")

  cat(paste("Downloading", data, "data from EBI AMR portal (ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/)\n"))

  if (!is.null(release)) {
    user_release <- release
    cat(paste(" Requesting data from user-specified release", user_release, "\n"))
  } else {
    # get list of releases
    folders <- str_split(RCurl::getURL(ebi_url, dirlistonly = TRUE), "\n")[[1]]
    # get latest
    user_release <- folders[!folders %in% c("releases.yml", "")] %>% max()
    cat(paste("...Requesting data from latest release", user_release, "\n"))
  }

  ebi_dat <- RCurl::getBinaryURL(print(gsub("ftp\\:", "https\\:", paste0(ebi_url, user_release, "/", paste0(data, ".parquet")))))

  cat(paste("...Reading parquet file\n"))
  ebi_dat <- arrow::read_parquet(ebi_dat)

  # filter by genus or species as required
  if (!is.null(genus)) {
    user_genus <- genus
    cat(paste("...Filtering by genus", genus, "\n"))
    ebi_dat <- ebi_dat %>%
      filter(genus == user_genus)
  } else if (!is.null(species)) {
    cat(paste("...Filtering by species", species, "\n"))
    user_species <- species
    ebi_dat <- ebi_dat %>%
      filter(species == user_species)
  }

  # filter by subclass/class/antibiotic as required
  if (!is.null(geno_subclass)) {
    if (data == "genotype") {
      cat(paste("...Filtering by subclass", paste(geno_subclass, collapse = ","), "\n"))
      ebi_dat <- ebi_dat %>%
        filter(grepl(paste(geno_subclass, collapse = "|"), subclass))
      for (sc in geno_subclass) {
        if (sum(grepl(sc, ebi_dat$subclass)) == 0) {
          cat(paste("......subclass not found:", sc, "\n"))
        }
      }
      if (!is.null(geno_class)) {
        cat(paste("...Warning, not using classes provided, filtering on subclass instead\n"))
      }
      if (!is.null(antibiotic)) {
        cat(paste("...Warning, not using antibiotic names provided, filtering on subclass instead\n"))
      }
    } else {
      cat(paste("...Warning, not using subclasses provided as these are only relevant to genotype data\n"))
    }
  } else if (!is.null(geno_class)) {
    if (data == "genotype") {
      cat(paste("...Filtering by class", paste(geno_class, collapse = ","), "\n"))
      ebi_dat <- ebi_dat %>%
        filter(grepl(paste(geno_class, collapse = "|"), class))
      for (cl in geno_subclass) {
        if (sum(grepl(cl, ebi_dat$class)) == 0) {
          cat(paste("......class not found:", cl, "\n"))
        }
      }
      if (!is.null(antibiotic)) {
        cat(paste("...Warning, not using antibiotic names provided, filtering on class instead\n"))
      }
    } else {
      cat(paste("...Warning, not using classes provided as these are only relevant to genotype data\n"))
    }
  } else if (!is.null(antibiotic)) {
    if (!force_antibiotic) {
      antibiotic <- na.omit(tolower(AMR::ab_name(AMR::as.ab(antibiotic))))
    }
    cat(paste("...Filtering by antibiotic:", paste(antibiotic, collapse = ", "), "\n"))
    ebi_dat <- ebi_dat %>%
      filter(antibiotic_name %in% antibiotic)
    for (ab in antibiotic) {
      if (!(ab %in% ebi_dat$antibiotic_name)) {
        cat(paste("......antibiotic not found:", ab, "\n"))
      }
    }
  }

  if (remove_dup) {
    cat("...Checking for duplicate rows ")
    before <- nrow(ebi_dat)
    ebi_dat <- ebi_dat %>%
      select(BioSample_ID:subclass) %>%
      distinct()
    after <- nrow(ebi_dat)
    if (after < before) {
      n <- before - after
      cat(paste("(removed", n, "rows)\n"))
    } else {
      cat(paste("(none detected)\n"))
    }
  }

  # reformat as per AMRgen import functions
  if (reformat) {
    if (data == "phenotype") {
      cat(paste("Reformatting phenotype data for easy use with AMRgen functions\n"))
      ebi_dat <- import_ebi_ast_ftp(ebi_dat,
        interpret_eucast = interpret_eucast,
        interpret_clsi = interpret_clsi,
        interpret_ecoff = interpret_ecoff
      )
    } else if (data == "genotype") {
      cat(paste("Reformatting genotype data for easy use with AMRgen functions\n"))
      ebi_dat <- import_amrfp_ebi(ebi_dat)
    }
  }

  # return data frame
  return(ebi_dat)
}
