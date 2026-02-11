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
#' Download NCBI antimicrobial susceptibility testing (AST) data
#'
#' This function downloads antimicrobial susceptibility testing (AST) data
#' from the NCBI Pathogen Detection database via the BioSample API. Data are
#' retrieved in batches, parsed from XML, and returned as a tidy tibble with
#' metadata including BioSample ID, Bioproject ID, and organism name.
#'
#' @param species Character. Organism name for the search query
#'   (e.g., `"Salmonella enterica"`). Required.
#' @param antibiotic Character or vector. Optional antibiotic name/s to filter
#'   the returned data. Strings will be processed using the \pkg{AMR} package
#'   to standardize names before matching, so e.g. `"amikacin"` or `"Amikacin"`
#'   or `"ami"` will be parsed to "amikacin" before matching. This can be turned
#'   off by setting `force_antibiotic=TRUE`. Full list of allowed antibiotic
#'   names in NCBI: <https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/>.
#' @param max_records Integer. Maximum number of BioSample records to
#'   retrieve. Default is `15000`.
#' @param batch_size Integer. Number of records fetched per API request.
#'   Default is `200` which is recommended by NCBI.
#' @param sleep_time Numeric. Seconds to pause between batch requests to
#'   avoid overloading NCBI servers. Default is `0.34`.
#' @param force_antibiotic Logical. If `TRUE`, turns off standardizing the
#'   antibiotic name using the \pkg{AMR} package before filtering, so that
#'   matching is done exactly on the input string/s. Default is `FALSE`.
#' @param reformat Logical. If `TRUE`, reformats the output using
#'   [import_ncbi_biosample] for compatibility with AMR analysis workflows.
#'   Default is `FALSE`. When set to `TRUE`, the data can also be interpreted
#'   against breakpoints/ECOFF by setting the `interpret_*=TRUE`.
#' @param interpret_eucast Logical. Passed to [interpret_ast] via
#'   [import_ncbi_biosample].
#'   If `TRUE`, interprets MIC values using EUCAST breakpoints. Default is
#'   `FALSE`. Only used if `reformat`=`TRUE`.
#' @param interpret_clsi Logical. Passed to [interpret_ast] via
#'   [import_ncbi_biosample].
#'   If `TRUE`, interprets MIC values using CLSI breakpoints. Default is
#'   `FALSE`. Only used if `reformat`=`TRUE`.
#' @param interpret_ecoff Logical. Passed to [interpret_ast] via
#'   [import_ncbi_biosample].
#'   If `TRUE`, interprets MIC values using ECOFF cutoffs. Default is `FALSE`.
#'   Only used if `reformat`=`TRUE`.
#'
#' @details
#' The function constructs an Entrez query of the form:
#' `"<organism> AND antibiogram[filter]"`. XML records are downloaded
#' in batches, parsed, and combined into a single table. The resulting tibble
#' contains AST test results and associated metadata including:
#'
#' - `id`: BioSample identifier
#' - `BioProject`: BioProject accession ID
#' - `organism`: Organism name
#' - `Antibiotic`, `Phenotype`, `Measurement`, `Units`, `Method`, `System`,
#'   `Manufacturer`, `Panel`, `Standard`: AST data columns
#'
#' The function can optionally filter by a one or more antibiotics. It can
#' also optionally reformat data for compatibility with AMRgen functions via
#' [import_ncbi_ast], and interpret the raw data measures against breakpoints or
#' ECOFF. See [import_ncbi_ast] for details of output formats when these options
#' are used.
#'
#' @return
#' A tibble with one row per AST measure, with corresponding BioSample metadata.
#'
#' @section NCBI API usage:
#' Users are encouraged to set an NCBI API key via
#' `rentrez::set_entrez_key()` to increase request limits and
#' comply with NCBI usage policies.
#'
#' @examples
#' \dontrun{
#' # Download AST data for Klebsiella quasipneumoniae
#' ast <- download_ncbi_ast("Klebsiella quasipneumoniae")
#'
#' # Download Klebsiella quasipneumoniae data, filter to amikacin and ampicillin
#' ast <- download_ncbi_ast(
#'   "Klebsiella quasipneumoniae",
#'   antibiotic = c("amikacin", "Amp")
#' )
#'
#' # Download and reformat for AMRgen workflow with EUCAST interpretation
#' ast <- download_ncbi_ast(
#'   "Klebsiella quasipneumoniae",
#'   reformat = TRUE,
#'   interpret_eucast = TRUE
#' )
#' }
#'
#' @import rentrez
#' @import tibble
#' @import XML
#' @importFrom AMR ab_name as.ab
#' @importFrom purrr modify_tree pluck
#' @importFrom dplyr bind_rows select filter mutate
#' @export
download_ncbi_ast <- function(species,
                              antibiotic = NULL,
                              max_records = 15000,
                              batch_size = 200,
                              sleep_time = 0.34,
                              force_antibiotic = FALSE,
                              reformat = FALSE,
                              interpret_eucast = FALSE,
                              interpret_clsi = FALSE,
                              interpret_ecoff = FALSE) {
  if (is.null(species)) {
    stop("`species` must be provided.")
  }

  # Build query term for entrez
  entrez_term <- paste0(stringr::str_trim(tolower(species)), "[orgn] AND antibiogram[filter]")

  # Search for AST entries by pathogen
  search <- rentrez::entrez_search(
    db = "biosample",
    term = entrez_term,
    retmax = max_records
  )

  # list samples by internal id
  ids <- search$ids

  if (length(ids) == 0) {
    warning("No AST records found for: ", species)
    return(tibble())
  }

  cat(paste(
    "Identified", length(ids), " ", tolower(species),
    " records from NCBI (https://www.ncbi.nlm.nih.gov/pathogens/ast/) \n"
  ))

  # create empty data structures to store records
  all_ast_data <- NULL

  # iterate though samples and pull records
  for (sample in seq(1, length(ids), by = batch_size)) {
    # account for remainder when batching
    if ((sample - 1 + batch_size > length(ids))) {
      data_xml <- rentrez::entrez_fetch(
        db = "biosample",
        id = ids[sample:(length(ids))],
        rettype = "xml",
        parsed = TRUE
      )

      cat(paste(
        "Downloading and processing records ", sample, " to ",
        length(ids), "... \n"
      ))
    } else {
      data_xml <- rentrez::entrez_fetch(
        db = "biosample",
        id = ids[sample:(sample - 1 + batch_size)],
        rettype = "xml",
        parsed = TRUE
      )

      cat(paste(
        "Downloading and processing records ", sample, " to ",
        (sample - 1 + batch_size), "... \n"
      ))
    }

    # Flatten XML
    data_list <- XML::xmlToList(data_xml)

    # extract records
    for (record in 1:length(data_list)) {
      data_names <- as.character(unlist(
        data_list[record]$BioSample$Description$Comment$Table$Header
      ))

      data_entries <- data_list[record]$BioSample$Description$Comment$Table$Body

      # retrieve AST data and reformat
      for (entry in 1:length(data_entries)) {
        temp_entry <- data_list[record]$BioSample$Description$Comment$Table$Body[entry]

        # Preserve null data points by converting to NA (for headers)
        temp_entry <- temp_entry %>%
          purrr::modify_tree(
            leaf = \(x) if (is.null(x)) NA else x,
            post = unlist
          ) %>%
          t() %>%
          as.data.frame()

        # add col names here to prevent data mismatches when binding rows
        names(temp_entry) <- c(as.character(data_names))

        # Extract BioProject accession while accounting for two varying flattened
        # xml node structures with duplicated identifiers
        bioproj <- purrr::pluck(data_list[record], "BioSample", "Links", 2, "Link", ".attrs", 3, .default = NA)

        if (is.na(bioproj)) {
          bioproj <- purrr::pluck(data_list[record], "BioSample", "Links", "Link", ".attrs", 3, .default = NA)
        } else {
          bioproj <- NA # when no data listed - very rare
        }

        # Extract organism name while accounting for two varied flattened
        # xml node structures
        organism_name <- purrr::pluck(data_list[record], "BioSample", "Description", "Organism", "OrganismName", .default = NA)

        if (is.na(organism_name)) {
          organism_name <- purrr::pluck(data_list[record], "BioSample", "Description", "Organism", "taxonomy_name", .default = NA)
        }

        # Add accessions and organism cols to AST data
        temp_entry <- temp_entry %>%
          mutate(id = data_list[record]$BioSample$Ids$Id$text) %>%
          mutate(BioProject = bioproj) %>%
          mutate(organism = organism_name)

        # combine records
        all_ast_data <- bind_rows(all_ast_data, temp_entry)
      }
    }

    # pause between record pulls
    Sys.sleep(sleep_time)
    cat(paste("Pausing download for ", sleep_time, "seconds to avoid overburdening the NCBI server... \n"))
  }

  # name and rearrange cols
  cat(paste("Ordering data... \n"))

  all_ast_data <- all_ast_data %>%
    select(id, BioProject, organism, everything()) %>%
    as_tibble()

  # filter data by drug where specified by user
  if (!is.null(antibiotic)) {
    if (!force_antibiotic) {
      antibiotic <- na.omit(tolower(AMR::ab_name(AMR::as.ab(antibiotic))))
    }

    cat(paste("...Filtering by antibiotic:", paste(antibiotic, collapse = ", "), "\n"))

    all_ast_data <- all_ast_data %>%
      filter(Antibiotic %in% antibiotic)
  }

  # reformat as per AMRgen import functions
  if (reformat) {
    cat(paste("Reformatting phenotype data for easy use with AMRgen functions \n"))

    all_ast_data <- import_ncbi_biosample(all_ast_data,
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff
    )
  }

  # return data
  return(all_ast_data)
}
