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

#' Export NCBI BioSample Antibiogram
#'
#' Convert AMRgen long-format AST data to an
#' [NCBI BioSample Antibiogram](https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/)
#' submission file.
#'
#' @param data A data frame in AMRgen long format (e.g. output of
#'   [import_ast()] or [format_ast()]).
#'   Expected columns: `id`, `drug_agent`, and at least one phenotype
#'   column (see `pheno_col`). Optional columns: `mic`, `disk`,
#'   `method`, `guideline`, `platform`.
#' @param file File path for the output file (must end in `.txt` or `.tsv`).
#' @param overwrite Logical; overwrite an existing file? Default `FALSE`.
#' @param pheno_col Character string naming the column that contains
#'   SIR interpretations (class `sir`). Default `"pheno_provided"`.
#'
#' @details
#' When both `mic` and `disk` columns are present, MIC values are
#' preferred (more precise). Disk values are only used for rows where
#' MIC is `NA`.
#'
#' MIC strings (e.g. `"<=0.5"`, `">=32"`, `"4"`) are split into a
#' sign (`<=`, `>=`, `<`, `>`, or `=`) and a numeric value.
#'
#' Antibiotic names are converted to lowercase with combination
#' separators replaced by `"-"` (NCBI convention, e.g.
#' `"amoxicillin-clavulanic acid"`).
#'
#' @return The formatted data frame is returned invisibly. A
#'   tab-delimited UTF-8 text file is written to `file`.
#'
#' @importFrom AMR ab_name
#' @importFrom dplyr mutate if_else case_when select any_of
#' @importFrom stringr str_match
#' @export
#' @examples
#' \dontrun{
#' # Write out the ecoli_ast data to file in NCBI format
#' export_ncbi_biosample(ecoli_ast, "Ec_NCBI.tsv")
#'
#' # Download data from EBI, then write it out to file in NCBI format
#' ebi_kleb_quasipneumoniae <- download_ebi(species = "Klebsiella quasipneumoniae", reformat = T)
#' export_ncbi_biosample(ebi_kq, "Kq_NCBI.tsv")
#' }
export_ncbi_biosample <- function(data, file, overwrite = FALSE,
                                  pheno_col = "pheno_provided") {
  # --- input validation ---
  if (file.exists(file) && !overwrite) {
    stop("The file '", file, "' already exists and `overwrite` is set to `FALSE`.")
  }
  if (!grepl("[.](txt|tsv)$", file, ignore.case = TRUE)) {
    stop("`file` must have the file extension '.txt' or '.tsv'.")
  }

  required <- c("id", "drug_agent", pheno_col)
  missing_req <- setdiff(required, colnames(data))
  if (length(missing_req) > 0) {
    stop("Missing required column(s): ", paste(missing_req, collapse = ", "))
  }

  has_mic <- "mic" %in% colnames(data)
  has_disk <- "disk" %in% colnames(data)
  if (!has_mic && !has_disk) {
    warning("Neither 'mic' nor 'disk' column found; measurement fields will be empty.")
  }

  # --- build measurement columns ---
  if (has_mic) {
    mic_str <- as.character(data$mic)
    mic_sign <- stringr::str_match(mic_str, "^(<=?|>=?)")[, 2]
    mic_sign <- dplyr::if_else(is.na(mic_sign) & !is.na(mic_str) & mic_str != "NA", "=", mic_sign)
    mic_value <- stringr::str_match(mic_str, "([0-9./]+)$")[, 2]
  }

  if (has_disk) {
    disk_str <- as.character(data$disk)
    disk_sign <- dplyr::if_else(!is.na(data$disk), "=", NA_character_)
    disk_value <- disk_str
  }

  # Determine sign, value, units — prefer MIC when both present
  if (has_mic && has_disk) {
    m_sign <- dplyr::if_else(!is.na(data$mic), mic_sign, disk_sign)
    m_value <- dplyr::if_else(!is.na(data$mic), mic_value, disk_value)
    m_units <- dplyr::if_else(!is.na(data$mic), "mg/L", dplyr::if_else(!is.na(data$disk), "mm", NA_character_))
  } else if (has_mic) {
    m_sign <- mic_sign
    m_value <- mic_value
    m_units <- dplyr::if_else(!is.na(data$mic), "mg/L", NA_character_)
  } else if (has_disk) {
    m_sign <- disk_sign
    m_value <- disk_value
    m_units <- dplyr::if_else(!is.na(data$disk), "mm", NA_character_)
  } else {
    m_sign <- rep(NA_character_, nrow(data))
    m_value <- rep(NA_character_, nrow(data))
    m_units <- rep(NA_character_, nrow(data))
  }

  # --- SIR → text ---
  sir_vals <- as.character(data[[pheno_col]])
  resistance_phenotype <- dplyr::case_when(
    sir_vals == "S" ~ "susceptible",
    sir_vals == "I" ~ "intermediate",
    sir_vals == "R" ~ "resistant",
    TRUE ~ NA_character_
  )

  # --- antibiotic name (lowercase, "-" for combos) ---
  antibiotic <- tryCatch(
    gsub("/", "-", AMR::ab_name(data$drug_agent, tolower = TRUE), fixed = TRUE),
    error = function(e) {
      warning("Could not convert some drug_agent values to antibiotic names: ", e$message)
      as.character(data$drug_agent)
    }
  )

  na_ab <- is.na(antibiotic) & !is.na(data$drug_agent)
  if (any(na_ab)) {
    warning(
      "AMR::ab_name() returned NA for some drug_agent values: ",
      paste(unique(data$drug_agent[na_ab]), collapse = ", ")
    )
  }

  # --- assemble output ---
  out <- data.frame(
    sample_name = data$id,
    antibiotic = antibiotic,
    resistance_phenotype = resistance_phenotype,
    measurement_sign = m_sign,
    measurement = m_value,
    measurement_units = m_units,
    laboratory_typing_method = if ("method" %in% colnames(data)) {
      dplyr::if_else(is.na(data$method), "missing", as.character(data$method))
    } else {
      "missing"
    },
    testing_standard = if ("guideline" %in% colnames(data)) {
      as.character(data$guideline)
    } else {
      NA_character_
    },
    laboratory_typing_platform = if ("platform" %in% colnames(data)) {
      as.character(data$platform)
    } else {
      NA_character_
    },
    vendor = NA_character_,
    laboratory_typing_method_version_or_reagent = NA_character_,
    stringsAsFactors = FALSE
  )

  # --- write ---
  utils::write.table(
    x = out,
    file = file,
    append = FALSE,
    quote = TRUE,
    sep = "\t",
    na = "",
    row.names = FALSE,
    col.names = TRUE,
    fileEncoding = "UTF-8"
  )

  invisible(out)
}


#' Export EBI Antibiogram
#'
#' Convert AMRgen long-format AST data to an EBI antibiogram
#' submission file (see
#' [EBI COMPARE-AMR](https://github.com/EBI-COMMUNITY/compare-amr)).
#'
#' @param data A data frame in AMRgen long format (e.g. output of
#'   [import_ast()] or [format_ast()]).
#'   Expected columns: `id`, `drug_agent`, `spp_pheno`, and at least
#'   one phenotype column (see `pheno_col`). Optional columns: `mic`,
#'   `disk`, `method`, `guideline`, `platform`.
#' @param file File path for the output file.
#' @param overwrite Logical; overwrite an existing file? Default `FALSE`.
#' @param pheno_col Character string naming the column that contains
#'   SIR interpretations (class `sir`). Default `"pheno_provided"`.
#' @param sep Field separator for the output file. Default `"\t"`
#'   (tab-delimited). Use `","` for CSV.
#'
#' @details
#' Antibiotic names are in Title Case with `"/"` separating
#' combination agents (EBI convention, e.g.
#' `"Amoxicillin/clavulanic acid"`).
#'
#' Species names are derived from the `spp_pheno` column via
#' [AMR::mo_name()].
#'
#' @return The formatted data frame is returned invisibly. A file is
#'   written to `file`.
#'
#' @importFrom AMR ab_name mo_name
#' @importFrom dplyr if_else case_when
#' @importFrom stringr str_match
#' @export
export_ebi_antibiogram <- function(data, file, overwrite = FALSE,
                                   pheno_col = "pheno_provided",
                                   sep = "\t") {
  # --- input validation ---
  if (file.exists(file) && !overwrite) {
    stop("The file '", file, "' already exists and `overwrite` is set to `FALSE`.")
  }

  required <- c("id", "drug_agent", pheno_col, "spp_pheno")
  missing_req <- setdiff(required, colnames(data))
  if (length(missing_req) > 0) {
    stop("Missing required column(s): ", paste(missing_req, collapse = ", "))
  }

  has_mic <- "mic" %in% colnames(data)
  has_disk <- "disk" %in% colnames(data)
  if (!has_mic && !has_disk) {
    warning("Neither 'mic' nor 'disk' column found; measurement fields will be empty.")
  }

  # --- build measurement columns ---
  if (has_mic) {
    mic_str <- as.character(data$mic)
    mic_sign <- stringr::str_match(mic_str, "^(<=?|>=?)")[, 2]
    mic_sign <- dplyr::if_else(is.na(mic_sign) & !is.na(mic_str) & mic_str != "NA", "=", mic_sign)
    mic_value <- stringr::str_match(mic_str, "([0-9./]+)$")[, 2]
  }

  if (has_disk) {
    disk_str <- as.character(data$disk)
    disk_sign <- dplyr::if_else(!is.na(data$disk), "=", NA_character_)
    disk_value <- disk_str
  }

  # Determine sign, value, units — prefer MIC when both present
  if (has_mic && has_disk) {
    m_sign <- dplyr::if_else(!is.na(data$mic), mic_sign, disk_sign)
    m_value <- dplyr::if_else(!is.na(data$mic), mic_value, disk_value)
    m_units <- dplyr::if_else(!is.na(data$mic), "mg/L", dplyr::if_else(!is.na(data$disk), "mm", NA_character_))
  } else if (has_mic) {
    m_sign <- mic_sign
    m_value <- mic_value
    m_units <- dplyr::if_else(!is.na(data$mic), "mg/L", NA_character_)
  } else if (has_disk) {
    m_sign <- disk_sign
    m_value <- disk_value
    m_units <- dplyr::if_else(!is.na(data$disk), "mm", NA_character_)
  } else {
    m_sign <- rep(NA_character_, nrow(data))
    m_value <- rep(NA_character_, nrow(data))
    m_units <- rep(NA_character_, nrow(data))
  }

  # --- SIR → text ---
  sir_vals <- as.character(data[[pheno_col]])
  resistance_phenotype <- dplyr::case_when(
    sir_vals == "S" ~ "susceptible",
    sir_vals == "I" ~ "intermediate",
    sir_vals == "R" ~ "resistant",
    TRUE ~ NA_character_
  )

  # --- species name ---
  species <- tryCatch(
    AMR::mo_name(data$spp_pheno),
    error = function(e) {
      warning("Could not convert some spp_pheno values to species names: ", e$message)
      as.character(data$spp_pheno)
    }
  )

  # --- antibiotic name (Title Case, "/" for combos) ---
  antibiotic_name <- tryCatch(
    AMR::ab_name(data$drug_agent),
    error = function(e) {
      warning("Could not convert some drug_agent values to antibiotic names: ", e$message)
      as.character(data$drug_agent)
    }
  )

  na_ab <- is.na(antibiotic_name) & !is.na(data$drug_agent)
  if (any(na_ab)) {
    warning(
      "AMR::ab_name() returned NA for some drug_agent values: ",
      paste(unique(data$drug_agent[na_ab]), collapse = ", ")
    )
  }

  # --- assemble output ---
  out <- data.frame(
    biosample_id = data$id,
    species = species,
    antibiotic_name = antibiotic_name,
    ast_standard = if ("guideline" %in% colnames(data)) {
      as.character(data$guideline)
    } else {
      NA_character_
    },
    breakpoint_version = NA_character_,
    laboratory_typing_method = if ("method" %in% colnames(data)) {
      as.character(data$method)
    } else {
      NA_character_
    },
    measurement = m_value,
    measurement_units = m_units,
    measurement_sign = m_sign,
    resistance_phenotype = resistance_phenotype,
    platform = if ("platform" %in% colnames(data)) {
      as.character(data$platform)
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )

  # --- write ---
  utils::write.table(
    x = out,
    file = file,
    append = FALSE,
    quote = TRUE,
    sep = sep,
    na = "",
    row.names = FALSE,
    col.names = TRUE,
    fileEncoding = "UTF-8"
  )

  invisible(out)
}


#' Export AST Data
#'
#' Generic dispatcher that exports AMRgen long-format AST data to a
#' submission-ready file. Currently supports NCBI BioSample
#' Antibiogram and EBI Antibiogram formats.
#'
#' @param data A data frame in AMRgen long format (e.g. output of
#'   [import_ast()] or [format_ast()]).
#'   Expected columns: `id`, `drug_agent`, `spp_pheno`, and at least
#'   one phenotype column (see `pheno_col`). Optional columns: `mic`,
#'   `disk`, `method`, `guideline`, `platform`.
#' @param file File path for the output file.
#' @param format Target format: `"ncbi"` (default) or `"ebi"`.
#' @param overwrite Logical; overwrite an existing file? Default `FALSE`.
#' @param pheno_col Character string naming the column that contains
#'   SIR interpretations. Default `"pheno_provided"`.
#' @param ... Additional arguments passed to the format-specific
#'   export function.
#'
#' @return The formatted data frame is returned invisibly.
#'
#' @export
export_ast <- function(data, file, format = "ncbi", overwrite = FALSE,
                       pheno_col = "pheno_provided", ...) {
  format <- tolower(format)
  switch(format,
    ncbi = export_ncbi_biosample(data, file,
      overwrite = overwrite,
      pheno_col = pheno_col
    ),
    ebi = export_ebi_antibiogram(data, file,
      overwrite = overwrite,
      pheno_col = pheno_col, ...
    ),
    stop("Unsupported format '", format, "'. Use 'ncbi' or 'ebi'.")
  )
}
