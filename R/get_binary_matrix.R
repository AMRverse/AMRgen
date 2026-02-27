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

#' Get Binary Matrix of Genotype and Phenotype Data
#'
#' This function generates a binary matrix representing the resistance (R vs S/I) and nonwildtype (NWT vs WT, or R/I vs S) status for a given antibiotic, and presence or absence of genetic markers related to one or more specified drug classes. It takes as input separate tables for genotype and phenotype data, matches these according to a common identifier (either specified by column names or assuming the first column contains the ID), and filters the data according to the specified antibiotic and drug class criteria before creating a binary matrix. Suitable input files can be generated using [import_ast()] to import phenotype data, and [import_amrfp()] to import genotype data from AMRFinderPlus.
#' @param geno_table A data frame containing genotype data, in long form with one row per sample and genetic marker. Expected format is that output by [import_amrfp()] and must include a column labeled `drug_class` (indicating the antibiotic class associated with each marker), in addition to a column indicating the marker (column name specified via `marker_col`) and a column for sample identifiers (specified via `geno_sample_col`, otherwise it is assumed the first column contains identifiers).
#' @param pheno_table A data frame containing phenotype data, in long form with one row per sample, drug and assay result. Expected format is that output by [import_ast()] and must include a column `drug_agent` (indicating the drug agent, interpretable as AMR pkg class `ab`), in addition to a column for sample identifiers (specified via `pheno_sample_col`, otherwise it is assumed the first column contains identifiers), a column with the resistance interpretation (S/I/R, specified via `sir_col`), and optionally a column with the ECOFF interpretation (WT/NWT or S/R, specified via `ecoff_col`).
#' @param antibiotic A character string specifying the antibiotic of interest to filter phenotype data. The value must match one of the entries in the `drug_agent` column of `pheno_table` or be coercible to a match using [as.ab].
#' @param drug_class_list A character vector of drug classes to filter genotype data for markers related to the specified antibiotic. Markers in `geno_table` will be filtered based on whether their `drug_class` matches any value in this list.
#' @param geno_sample_col A character string (optional) specifying the column name in `geno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers.
#' @param pheno_sample_col A character string (optional) specifying the column name in `pheno_table` containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column contains identifiers.
#' @param sir_col A character string specifying the column name in `pheno_table` that contains the resistance interpretation (SIR) data. The values should be `"S"`, `"I"`, `"R"` or otherwise interpretable by [AMR::as.sir()]. If not provided, the first column prefixed with "phenotype*" will be used if present, otherwise an error is thrown.  Only used if `binary_matrix` not provided.
#' @param ecoff_col A character string specifying the column name in `pheno_table` that contains the ECOFF interpretation of phenotype. The values should be interpretable as `"WT"` (wildtype) and `"NWT"` (nonwildtype), or `"S"` / `"I"` / `"R"`. Default `"ecoff`.
#' @param marker_col A character string specifying the column name in `geno_table` containing the marker identifiers. Defaults to `"marker"`.
#' @param keep_SIR A logical indicating whether to retain the full S/I/R phenotype column in the output. Defaults to `TRUE`.
#' @param keep_assay_values A logical indicating whether to include columns with the raw phenotype assay data. Assumes there are columns labelled `"mic"` and `"disk"`; these will be added to the output table if present. Defaults to `FALSE`.
#' @param keep_assay_values_from A character vector specifying which assay values (e.g., `"mic"`, `"disk"`) to retain if `keep_assay_values` is `TRUE`. Defaults to `c("mic", "disk")`.
#' @param most_resistant A logical indicating whether, when multiple phenotype entries are present for the same sample and drug, the most resistant should be kept (otherwise the least resistant is kept). Default is `TRUE`.
#' @importFrom AMR as.ab as.disk as.mic as.sir is.ab
#' @importFrom dplyr any_of arrange case_when count desc filter full_join group_by if_else mutate pull right_join select slice_head ungroup
#' @importFrom tidyr pivot_wider
#' @return A data frame where each row represents a sample, and each column represents a genetic marker related to the specified antibiotic's drug class. The binary values in the matrix indicate the presence (`1`) or absence (`0`) of each marker for each sample, along with resistance status columns for the specified antibiotic: `R` for resistant (defined from `sir_col`, 1=R, 0=I/S) and `NWT` for nonwildtype (defined by `ecoff_col` if provided: 1=NWT, 0=WT; otherwise defined from `sir_col`: 1=I/R, 0=S).
#' @details This function performs several steps:
#' - Verifies that the `pheno_table` contains a `drug_agent` column and converts it to class `ab` if necessary.
#' - Filters the `pheno_table` to retain data related to the specified antibiotic.
#' - Checks that the `geno_table` contains markers associated with the specified drug class(es).
#' - Matches sample identifiers between `geno_table` and `pheno_table`.
#' - Extracts and transforms the phenotype data into a binary format indicating resistance and NWT status.
#' - Constructs a binary matrix for genotype data, with each column representing a genetic marker.
#' - Returns a single matrix with sample identifiers plus binary variables for each phenotype and genetic marker.
#' @seealso `compare_geno_pheno_id`, `as.ab`, `as.sir`, `ab_name`
#' @export
#' @examples
#' \dontrun{
#' # Import example E. coli AMRFinderPlus data from AllTheBacteria
#' ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
#' geno_pheno_cip <- get_binary_matrix(
#'   ecoli_geno,
#'   ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno_clsi"
#' )
#' geno_pheno_cip <- get_binary_matrix(
#'   ecoli_geno,
#'   ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "Resistance phenotype",
#'   keep_assay_values = TRUE
#' )
#' }
get_binary_matrix <- function(geno_table, pheno_table, antibiotic, drug_class_list, keep_SIR = TRUE,
                              keep_assay_values = FALSE, keep_assay_values_from = c("mic", "disk"),
                              geno_sample_col = NULL, pheno_sample_col = NULL,
                              sir_col = "pheno_clsi", ecoff_col = "ecoff", marker_col = "marker",
                              most_resistant = TRUE) {
  
  # check there is a SIR column specified
  if (is.null(sir_col)) {
    # make a sensible guess
    sir_col <- pheno_table %>%
      select(starts_with("pheno")) %>%
      colnames() %>%
      first()
    if (!is.na(sir_col)) {
      cat(paste("WARNING: `sir_col` not provided, using first column with prefix 'pheno':", sir_col, "\n"))
    } else {
      stop("`sir_col` not provided. Please specify a column with S/I/R phenotype values.")
    }
  }
  if (!(sir_col %in% colnames(pheno_table))) {
    stop(paste0("Column: '", sir_col, "' not found in input phenotype data. Please specify a valid column with S/I/R values."))
  }
  
  # check we have a drug_agent column with class ab
  if (!("drug_agent" %in% colnames(pheno_table))) {
    stop(paste("input", deparse(substitute(pheno_table)), "must have a column labelled `drug_agent`"))
  }
  if (!is.ab(pheno_table$drug_agent)) {
    cat(paste(" Converting", deparse(substitute(pheno_table)), "column `drug_agent` to class `ab` in binary matrix\n"))
    pheno_table <- pheno_table %>% mutate(drug_agent = as.ab(drug_agent))
  }

  # filter pheno to the drug of interest
  if (as.ab(antibiotic) %in% pheno_table$drug_agent) {
    pheno_table <- pheno_table %>% filter(drug_agent == as.ab(antibiotic))
  } else {
    stop(paste("antibiotic ", antibiotic, " was not found in input: ", deparse(substitute(pheno_table))))
  }

  # check we have some geno hits for markers relevant to the drug class/es
  if (!("drug_class" %in% colnames(geno_table))) {
    stop(paste("input", deparse(substitute(geno_table)), "must have a column labelled `drug_class`"))
  }
  if (sum(drug_class_list %in% geno_table$drug_class) == 0) {
    stop(paste("No markers matching drug class", paste(drug_class_list, collapse = ","), "were identified in input geno_table"))
  }

  # subset pheno & geno dataframes to those samples with overlap
  overlap <- compare_geno_pheno_id(geno_table, pheno_table, geno_sample_col = geno_sample_col, pheno_sample_col = pheno_sample_col, rename_id_cols = TRUE)
  pheno_matched <- overlap$pheno_matched
  geno_matched <- overlap$geno_matched

  # take single representative phenotype row per sample
  pheno_matched_rows_unfiltered <- nrow(pheno_matched)

  if (most_resistant) { # take most resistant
    pheno_matched <- pheno_matched %>%
      arrange(desc(get(sir_col)), desc(mic)) %>%
      group_by(id) %>%
      slice_head(n = 1) %>%
      ungroup()
    if (nrow(pheno_matched) < pheno_matched_rows_unfiltered) {
      cat(" Some samples had multiple phenotype rows, taking the most resistant only for binary matrix\n")
    }
  }
  if (!most_resistant) { # take least resistant
    pheno_matched_rows_unfiltered <- nrow(pheno_matched)
    pheno_matched <- pheno_matched %>%
      arrange(get(sir_col), mic) %>%
      group_by(id) %>%
      slice_head(n = 1) %>%
      ungroup()
    if (nrow(pheno_matched) < pheno_matched_rows_unfiltered) {
      cat("Some samples had multiple phenotype rows, taking the least resistant only\n")
    }
  }

  # check we have retained some samples that have relevant phenotype, and genotype data
  if (nrow(pheno_matched) == 0 | nrow(geno_matched) == 0) {
    stop(paste("No samples with both phenotype data for", antibiotic, "and genotype results"))
  }

  # get interpreted phenotype as binary (based on colname provided by 'sir_col')
  # to do: check we have interpretation data for this pheno, and optionally interpret from mic/disk
  pheno_binary <- pheno_matched %>%
    select(id, any_of(c(sir_col, ecoff_col))) %>%
    mutate(R = case_when(
      as.sir(get(sir_col)) == "R" ~ 1,
      as.sir(get(sir_col)) == "I" ~ 0,
      as.sir(get(sir_col)) == "S" ~ 0,
      TRUE ~ NA
    )) %>%
    mutate(NWT = case_when(
      as.sir(get(sir_col)) == "R" ~ 1,
      as.sir(get(sir_col)) == "I" ~ 1,
      as.sir(get(sir_col)) == "S" ~ 0,
      TRUE ~ NA
    ))

  # replace NWT with ecoff-based definition if available
  if (!is.null(ecoff_col)) {
    if (ecoff_col %in% colnames(pheno_binary)) {
      cat(paste(" Defining NWT in binary matrix using ecoff column provided:", ecoff_col, "\n"))
      pheno_binary <- pheno_binary %>%
        mutate(NWT = case_when(
          as.sir(get(ecoff_col)) == "NWT" ~ 1,
          as.sir(get(ecoff_col)) == "WT" ~ 0,
          as.sir(get(ecoff_col)) == "R" ~ 1,
          as.sir(get(ecoff_col)) == "I" ~ 1,
          as.sir(get(ecoff_col)) == "S" ~ 0,
          TRUE ~ NA
        ))
    } else {
      cat(" Defining NWT in binary matrix as I/R vs 0, as no ECOFF column defined\n")
    }
  }

  pheno_binary <- pheno_binary %>% select(id, R, NWT)

  # check there are some non-NA values for R/NWT calls
  if (sum(!is.na(pheno_binary$R)) == 0 & sum(!is.na(pheno_binary$NWT)) == 0) {
    stop(paste("No samples with both genotype data and non-NA phenotype interpretation values for", antibiotic, "in input", deparse(substitute(pheno_table))))
  }

  # extract list of relevant drug markers
  markers <- geno_matched %>%
    filter(drug_class %in% drug_class_list | as.ab(drug_agent) == as.ab(antibiotic)) %>%
    pull(get(marker_col)) %>%
    unique()

  # get geno as binary table indicating presence/absence for the relevant drug markers
  geno_binary <- geno_matched %>%
    filter(get(marker_col) %in% markers) %>%
    group_by(id, get(marker_col)) %>%
    count() %>%
    rename(marker = `get(marker_col)`) %>%
    mutate(n = if_else(n > 1, 1, n)) %>% # only count 1 per strain
    ungroup() %>%
    right_join(pheno_binary, by = "id") %>%
    mutate(marker = gsub(":", "..", marker)) %>% # replace : separator with .. otherwise pivot_wider replaces it with _
    pivot_wider(names_from = marker, values_from = n, values_fill = 0)

  # if there were samples with phenotypes, but no hits for any markers, there will be a 'NA' column created, need to remove this
  if ("NA" %in% colnames(geno_binary)) {
    geno_binary <- geno_binary %>% select(-c("NA"))
  }

  # add MIC / disk values if specified
  if (keep_assay_values) {
    if (sum(keep_assay_values_from %in% colnames(pheno_matched)) > 0) {
      geno_binary <- pheno_matched %>%
        select(id, any_of(keep_assay_values_from)) %>%
        full_join(geno_binary, by = "id")
    } else {
      cat(paste(" No specified assay columns found to include in binary matrix:", keep_assay_values_from, "\n"))
    }
    if ("mic" %in% colnames(geno_binary)) {
      geno_binary <- geno_binary %>% mutate(mic = as.mic(mic))
    }
    if ("disk" %in% colnames(geno_binary)) {
      geno_binary <- geno_binary %>% mutate(disk = as.disk(disk))
    }
  }

  # add SIR values unless switched off
  if (keep_SIR) {
    sir_binary <- pheno_matched %>%
      select(id, any_of(c(sir_col, ecoff_col))) %>%
      mutate(pheno = as.sir(get(sir_col)))

    if (!is.null(ecoff_col)) {
      if (ecoff_col %in% colnames(sir_binary)) {
        sir_binary <- sir_binary %>% mutate(ecoff = as.sir(get(ecoff_col)))
      }
    }

    geno_binary <- sir_binary %>%
      select(id, any_of(c("pheno", "ecoff"))) %>%
      full_join(geno_binary, by = "id")
  }

  return(geno_binary)
}


#' Add marker combinations to a binary geno-pheno matrix
#'
#' Given a geno-pheno binary marker matrix, output by `get_binary_matrix`,
#' this function identifies marker combination, and adds combination ids,
#' and the count of markers detected per sample, as new columns.
#'
#' @param binary_matrix A geno-pheno binary matrix, output by `get_binary_matrix`
#' @param assay (optional) Name of an assay column to filter on, so that the matrix returned only includes samples with assay data of this type available
#'
#'
#' @examples
#' \dontrun{
#' combo_matrix <- get_combo_matrix(binary_matrix)
#' }
#'
#' @export
get_combo_matrix <- function(binary_matrix, assay = NULL) {
  # marker names
  markers <- binary_matrix %>%
    dplyr::select(-any_of(c("id", "pheno", "ecoff", "R", "I", "NWT", "mic", "disk"))) %>%
    colnames()

  # check we have the expected assay column, with data
  if (!is.null(assay)) {
    if (!(assay %in% colnames(binary_matrix))) {
      stop(paste("input", deparse(substitute(binary_matrix)), "has no column labelled ", assay))
    }
    data_rows <- binary_matrix %>%
      filter(!is.na(get(assay))) %>%
      nrow()
    if (data_rows == 0) {
      stop(paste("input", deparse(substitute(binary_matrix)), "has no non-NA values in column ", assay))
    }
  }

  # Add in a combination column and filter to samples with the required assay data (MIC or disk) only
  if (!is.null(assay)) {
    binary_matrix <- binary_matrix %>% filter(!is.na(get(assay)))
  }

  # input matrix with combinations added
  binary_matrix_combo <- binary_matrix %>%
    mutate(marker_count = rowSums(across(where(is.numeric) & !any_of(c("R", "NWT", "mic", "disk"))), na.rm = T)) %>%
    unite("combination_id", markers[1]:markers[length(markers)], remove = FALSE) %>% # add in combinations
    relocate(marker_count, .after = "combination_id")

  return(binary_matrix_combo)
}
