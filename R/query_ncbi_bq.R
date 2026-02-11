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

#' Query antimicrobial phenotype (AST) data from NCBI Pathogen Detection BigQuery
#'
#' This function queries the `ncbi-pathogen-detect.pdbrowser.ast` BigQuery table to retrieve
#' antimicrobial susceptibility testing (AST) data from NCBI Pathogen Detection.
#'
#' @details
#' Requires Google Cloud authentication. Run `bigrquery::bq_auth()` before
#' first use, or set up application default credentials via
#' `gcloud auth application-default login`.
#'
#' @param taxgroup String specifying the organism group to filter on (e.g.,
#' "Pseudomonas aeruginosa"). See https://www.ncbi.nlm.nih.gov/pathogens/organisms/
#' for a list. Required.
#' @param antibiotic (Optional) String (or vector of strings) specifying the
#' antibiotic name/s to filter on (default NULL). Uses the AMR package to try
#' to fix typos, and format to lower-case.
#' @param force_antibiotic (Optional) Logical indicating whether to turn off
#' parsing of antibiotic names and match exactly on the input strings (default
#' FALSE).

#' @param project_id (Optional) Google Cloud Project ID to use for billing. If
#'  NULL (default), looks for `GOOGLE_CLOUD_PROJECT` environment variable.
#'
#' @return A tibble containing AST data with columns renamed to match
#' `import_ncbi_ast()` expectations.
#' @importFrom bigrquery bq_project_query bq_table_download
#' @importFrom dplyr rename mutate
#' @export
#' @examples
#' \dontrun{
#' # Query AST data for Klebsiella pneumoniae, filtered to meropenem
#' ast_raw <- query_ncbi_bq_ast(
#'     taxgroup = "Klebsiella pneumoniae",
#'     antibiotic = "meropenem"
#' )
#'
#' # Import and reinterpret using CLSI breakpoints
#' ast <- import_ncbi_ast(ast_raw, interpret_clsi = TRUE)
#' }
query_ncbi_bq_ast <- function(taxgroup,
                              antibiotic = NULL,
                              force_antibiotic = FALSE,
                              project_id = NULL) {
    if (missing(taxgroup)) stop("Argument 'taxgroup' is required.")

    project_id <- get_bq_project_id(project_id)

    # Base query
    query <- "
    SELECT
      biosample_acc,
      platform,
      standard,
      disk_diffusion,
      mic,
      measurement_sign,
      antibiotic,
      scientific_name,
      bioproject_acc,
      phenotype
    FROM `ncbi-pathogen-detect.pdbrowser.ast`
    WHERE taxgroup_name = @taxgroup
  "

    params <- list(taxgroup = taxgroup)


    # Add antibiotic filter
    if (!is.null(antibiotic)) {
        if (!force_antibiotic) {
            if (requireNamespace("AMR", quietly = TRUE)) {
                antibiotic <- unique(tolower(AMR::ab_name(AMR::as.ab(antibiotic))))
            } else {
                warning("AMR package not installed. Using provided antibiotic names without parsing.")
                antibiotic <- tolower(antibiotic)
            }
        } else {
            antibiotic <- tolower(antibiotic)
        }

        if (length(antibiotic) == 1) {
            query <- paste0(query, " AND antibiotic = @antibiotic")
        } else {
            query <- paste0(query, " AND antibiotic IN UNNEST(@antibiotic)")
        }
        params$antibiotic <- antibiotic
    }

    # Execute query
    res <- bq_query_with_auth_check(project_id, query, params)

    # Rename columns to match import_ncbi_ast expectations
    res <- res %>%
        dplyr::rename(
            "BioSample" = biosample_acc,
            "Laboratory typing platform" = platform,
            "Testing standard" = standard,
            "Disk diffusion (mm)" = disk_diffusion,
            "MIC (mg/L)" = mic,
            "Measurement sign" = measurement_sign,
            "Antibiotic" = antibiotic,
            "BioProject" = bioproject_acc,
            "Resistance phenotype" = phenotype,
            "Scientific name" = scientific_name
        )

    return(res)
}

#' Query antimicrobial genotype (MicroBIGG-E) data from NCBI Pathogen Detection BigQuery
#'
#' This function queries the `ncbi-pathogen-detect.pdbrowser.microbigge` BigQuery table to retrieve
#' genotype data. **Note:** This function only returns genotypes for BioSamples that also have AST data.
#'
#' @inheritParams query_ncbi_bq_ast
#' @param geno_subclass (Optional) String or vector of strings specifying AMR
#' subclasses to filter on (e.g., "CARBAPENEM").
#' @param geno_class (Optional) String or vector of strings specifying AMR
#' classes to filter on (e.g., "BETA-LACTAM"). Ignored if `geno_subclass` is
#' specified.
#'
#' @return A tibble containing genotype data with columns renamed to match
#' `import_amrfp()` expectations.
#' @importFrom bigrquery bq_project_query bq_table_download
#' @importFrom dplyr rename
#' @export
#' @examples
#' \dontrun{
#' # Query genotype data for Klebsiella pneumoniae samples that have AST data
#' geno_raw <- query_ncbi_bq_geno(
#'     taxgroup = "Klebsiella pneumoniae"
#' )
#'
#' # Filter for carbapenem resistance genes and point mutations
#' geno_amrfp <- query_ncbi_bq_geno(
#'     taxgroup = "Klebsiella pneumoniae",
#'     geno_subclass = "CARBAPENEM"
#' )
#'
#' geno <- import_amrfp(geno_amrfp, sample_col = "biosample_acc")
#' }
query_ncbi_bq_geno <- function(taxgroup,
                               geno_subclass = NULL,
                               geno_class = NULL,
                               project_id = NULL) {
    if (missing(taxgroup)) stop("Argument 'taxgroup' is required.")

    project_id <- get_bq_project_id(project_id)

    # Base query with CTE to filter for samples with AST data
    query <- "
    WITH has_ast as (
      SELECT distinct(biosample_acc)
      FROM `ncbi-pathogen-detect.pdbrowser.ast`
      WHERE taxgroup_name = @taxgroup
  "
    params <- list(taxgroup = taxgroup)

    query <- paste0(query, ") ")

    # Main query
    query <- paste0(query, "
    SELECT
      mb.biosample_acc,
      element_symbol,
      subclass,
      type,
      subtype,
      amr_method,
      hierarchy_node,
      scientific_name
    FROM
    `ncbi-pathogen-detect.pdbrowser.microbigge` mb
    INNER JOIN has_ast on has_ast.biosample_acc = mb.biosample_acc
    WHERE type = 'AMR'
  ")

    # Add subclass/class filters
    if (!is.null(geno_subclass)) {
        if (length(geno_subclass) == 1) {
            query <- paste0(query, " AND subclass = @subclass")
        } else {
            query <- paste0(query, " AND subclass IN UNNEST(@subclass)")
        }
        params$subclass <- geno_subclass
    } else if (!is.null(geno_class)) {
        if (length(geno_class) == 1) {
            query <- paste0(query, " AND class = @class")
        } else {
            query <- paste0(query, " AND class IN UNNEST(@class)")
        }
        params$class <- geno_class
    }

    # Execute query
    res <- bq_query_with_auth_check(project_id, query, params)

    # Rename columns to match import_amrfp expectations
    res <- res %>%
        dplyr::rename(
            "Gene symbol" = element_symbol,
            "Element type" = type,
            "Method" = amr_method,
            "Hierarchy_node" = hierarchy_node,
            "Element subtype" = subtype,
            "Subclass" = subclass
        )

    return(res)
}

# Internal helper to resolve project ID
get_bq_project_id <- function(project_id = NULL) {
    if (!is.null(project_id)) {
        return(project_id)
    }

    env_project <- Sys.getenv("GOOGLE_CLOUD_PROJECT", unset = "")
    if (nzchar(env_project)) {
        return(env_project)
    }

    stop("GCP Project ID must be provided via `project_id` argument or `GOOGLE_CLOUD_PROJECT` environment variable. \nSee https://cloud.google.com/resource-manager/docs/creating-managing-projects#identifying_projects", call. = FALSE)
}

# Internal helper to execute query with auth error handling
bq_query_with_auth_check <- function(project_id, query, params) {
    tryCatch(
        {
            job <- bigrquery::bq_project_query(project_id, query = query, parameters = params)
            bigrquery::bq_table_download(job)
        },
        error = function(e) {
            if (grepl("authentication|credentials|401", e$message, ignore.case = TRUE)) {
                stop("BigQuery authentication failed. Please run `bigrquery::bq_auth()` to authenticate, or ensure application default credentials are set.", call. = FALSE)
            }
            stop(e)
        }
    )
}
