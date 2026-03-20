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
#'   Expected columns: `id`, `drug`, and at least one phenotype
#'   column (see `pheno_col`). Optional columns: `mic`, `disk`,
#'   `method`, `guideline`, `platform`.
#' @param file File path for the output file (must end in `.txt` or
#'   `.tsv`). If `NULL` (default), no file is written and the
#'   formatted data frame is returned visibly.
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
#' @return When `file` is provided, the formatted data frame is
#'   returned invisibly and a tab-delimited UTF-8 file is written to
#'   `file`. When `file = NULL`, the formatted data frame is returned
#'   visibly and no file is written.
#'
#' @importFrom AMR ab_name
#' @importFrom dplyr mutate if_else case_when select any_of
#' @importFrom stringr str_match
#' @importFrom readr write_tsv
#' @export
#' @examples
#' \dontrun{
#' # Return formatted data frame without writing a file
#' ncbi_df <- export_ncbi_ast(ecoli_ast)
#'
#' # Write out the ecoli_ast data to file in NCBI format
#' export_ncbi_ast(ecoli_ast, "Ec_NCBI.tsv")
#'
#' # Download data from EBI, then write it out to file in NCBI format
#' ebi_kq <- download_ebi(
#'   data = "phenotype",
#'   species = "Klebsiella quasipneumoniae",
#'   reformat = T
#' )
#' export_ncbi_ast(ebi_kq, "Kq_NCBI.tsv")
#' }
export_ncbi_ast <- function(data, file = NULL, overwrite = FALSE,
                            pheno_col = "pheno_provided") {
  # --- input validation ---
  if (!is.null(file)) {
    if (file.exists(file) && !overwrite) {
      stop("The file '", file, "' already exists and `overwrite` is set to `FALSE`.")
    }
    if (!grepl("[.](txt|tsv)$", file, ignore.case = TRUE)) {
      stop("`file` must have the file extension '.txt' or '.tsv'.")
    }
  }

  required <- c("id", "drug", pheno_col)
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
    gsub("/", "-", AMR::ab_name(data$drug, tolower = TRUE), fixed = TRUE),
    error = function(e) {
      warning("Could not convert some drug values to antibiotic names: ", e$message)
      as.character(data$drug)
    }
  )

  na_ab <- is.na(antibiotic) & !is.na(data$drug)
  if (any(na_ab)) {
    warning(
      "AMR::ab_name() returned NA for some drug values: ",
      paste(unique(data$drug[na_ab]), collapse = ", ")
    )
  }

  # --- map method to NCBI-permitted values ---
  # Permitted: "MIC", "agar dilution", "disk diffusion", "missing"
  method_raw <- if ("method" %in% colnames(data)) as.character(data$method) else rep(NA_character_, nrow(data))
  ncbi_method <- dplyr::case_when(
    method_raw %in% c("MIC", "broth dilution") ~ "MIC",
    method_raw == "disk diffusion" ~ "disk diffusion",
    method_raw == "agar dilution" ~ "agar dilution",
    TRUE ~ "missing"
  )

  # --- assemble output ---
  out <- data.frame(
    sample_name = data$id,
    antibiotic = antibiotic,
    resistance_phenotype = resistance_phenotype,
    measurement_sign = m_sign,
    measurement = m_value,
    measurement_units = m_units,
    laboratory_typing_method = ncbi_method,
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
  if (!is.null(file)) {
    write_tsv(
      x = out,
      file = file,
      append = FALSE,
      na = "",
      col_names = TRUE,
      quote = "all",
      eol = "\n"
    )
    return(invisible(out))
  }

  out
}


#' Export EBI Antibiogram
#'
#' Format AMRgen long-format AST data to a table with the fields required
#' for submission to EBI, and optionally generate JSON submission files (one
#' per BioSample).
#' See
#' [https://www.ebi.ac.uk/amr/amr_submission_guide/](https://www.ebi.ac.uk/amr/amr_submission_guide/)).
#'
#' @param data A data frame in AMRgen long format (e.g. output of
#'   [import_ast()] or [format_ast()]).
#'   Expected columns: `id`, `drug`, `spp_pheno`, and at least
#'   one phenotype column (see `pheno_col`). Optional columns: `mic`,
#'   `disk`, `method`, `platform`.
#' @param pheno_col Character string naming the column that contains
#'   SIR interpretations (class `sir`). Default `"pheno_provided"`.
#' @param breakpoint_version Character string specifying the breakpoint
#' version used for interpretation (e.g. `"EUCAST 2024"`).
#' @param submission_account Character string specifying the EBI Webin
#' submission account identifier (e.g. `"Webin-###"`). If not provided,
#' JSON output files will not be generated and the function will return
#' the formated table only, which can be further updated and converted to
#' submission-ready JSON later using [format_ebi_json()].
#' @param domain Character string specifying the domain used in the
#' submission metadata (e.g. `"self.ExampleDomain"`). If not provided,
#' JSON output files will not be generated and the function will return
#' the formated table only, which can be further updated and converted to
#' submission-ready JSON later using [format_ebi_json()].
#' @param output_dir Character string specifying the directory where JSON
#' files should be written. If not provided,
#' JSON output files will not be generated and the function will return
#' the formated table only, which can be further updated and converted to
#' submission-ready JSON later using [format_ebi_json()].
#'
#' @details
#' Antibiotic names are in Title Case with `"/"` separating
#' combination drugs (EBI convention, e.g.
#' `"Amoxicillin/clavulanic acid"`).
#'
#' Species names are derived from the `spp_pheno` column via
#' [AMR::mo_name()].
#'
#' @return Formatted data frame. When `output_dir` is provided, the AST data
#' is also written to individual JSON submission files, one per BioSample, in
#' the specified directory.
#'
#' @importFrom AMR ab_name mo_name
#' @importFrom dplyr if_else case_when
#' @importFrom stringr str_match
#' @importFrom readr write_tsv
#' @export
#' @examples
#' # Return formatted data frame without writing files
#' ebi_df <- export_ebi_ast(staph_ast_ebi)
#' \dontrun{
#' # Write out data for each BioSample to an individual JSON file for submission
#' ebi_df <- export_ebi_ast(staph_ast_ebi,
#'   breakpoint_version = "EUCAST 2015",
#'   submission_account = "Webin-###",
#'   domain = "self.ExampleDomain",
#'   output_dir = "/path/to/output/"
#' )
#' }
export_ebi_ast <- function(data,
                           pheno_col = "pheno_provided",
                           breakpoint_version,
                           submission_account,
                           domain = "self.ExampleDomain",
                           output_dir = NULL) {
  # --- input validation ---
  required <- c("id", "drug", pheno_col, "spp_pheno")
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
    AMR::ab_name(data$drug),
    error = function(e) {
      warning("Could not convert some drug values to antibiotic names: ", e$message)
      as.character(data$drug)
    }
  )

  na_ab <- is.na(antibiotic_name) & !is.na(data$drug)
  if (any(na_ab)) {
    warning(
      "AMR::ab_name() returned NA for some drug values: ",
      paste(unique(data$drug[na_ab]), collapse = ", ")
    )
  }

  # --- map method to EBI-permitted values ---
  # Permitted: "E-test", "agar dilution", "broth dilution", "disk diffusion"
  method_raw_ebi <- if ("method" %in% colnames(data)) as.character(data$method) else rep(NA_character_, nrow(data))
  ebi_method <- dplyr::case_when(
    method_raw_ebi %in% c("MIC", "broth dilution") ~ "broth dilution",
    method_raw_ebi == "disk diffusion" ~ "disk diffusion",
    method_raw_ebi == "agar dilution" ~ "agar dilution",
    tolower(method_raw_ebi) %in% c("etest", "e-test") ~ "E-test",
    TRUE ~ NA_character_
  )

  # --- assemble tabular output ---
  out_df <- data.frame(
    biosample_id = data$id,
    species = species,
    antibiotic_name = antibiotic_name,
    ast_standard = if ("guideline" %in% colnames(data)) {
      as.character(data$guideline)
    } else {
      NA_character_
    },
    breakpoint_version = NA_character_,
    laboratory_typing_method = ebi_method,
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
  if (!is.null(output_dir)) {
    safe_execute(format_ebi_json(out_df,
      output_dir = output_dir,
      breakpoint_version = breakpoint_version,
      submission_account = submission_account,
      domain = domain
    ))
    return(invisible(out_df))
  }

  out_df
}


#' Generate EBI antibiogram submission in JSON
#'
#' Converts the tabular output of [export_ebi_ast()] into JSON
#' files formatted for submission to EBI as BioSample data
#' (https://www.ebi.ac.uk/amr/amr_submission_guide/). Each row of the input
#' dataset is converted into JSON records and printed to file.
#'
#' @param ebi_antibiogram_table A data frame in the format output by
#' [export_ebi_ast()].
#' @param breakpoint_version Character string specifying the
#' breakpoint version used for interpretation (e.g. `"EUCAST 2024"`).
#' @param submission_account Character string specifying the Webin
#' submission account identifier (e.g. `"Webin-###"`).
#' @param output_dir Character string specifying the directory where JSON
#' files should be written.
#' @param domain (Optional) Character string specifying the domain used in the
#' submission metadata (default `"self.ExampleDomain"`).
#' @return Invisibly returns `NULL`. The function prints JSON-formatted
#' AMR submission records to file.
#'
#' @details
#' The function iterates over each biosample in `ebi_antibiogram_table` and
#' constructs
#' a nested JSON object describing the antimicrobial susceptibility
#' testing result. Each record contains antibiotic metadata, AST
#' standards, measurement values, and resistance phenotype information.
#'
#' JSON formatting is performed using \code{jsonlite::toJSON()} with
#' `pretty = TRUE` and `auto_unbox = TRUE` .
#'
#' @examples
#' \dontrun{
#' format_ebi_json(
#'   ast_dataset,
#'   breakpoint_version = "EUCAST 2015",
#'   submission_account = "Webin-###",
#'   output_dir = "/path/to/output/"
#' )
#' }
#'
#' @importFrom jsonlite write_json
#' @export
format_ebi_json <- function(ebi_antibiogram_table,
                            breakpoint_version,
                            submission_account,
                            output_dir,
                            domain = NULL) {
  if (!dir.exists(output_dir)) {
    safe_execute(dir.create(output_dir, recursive = TRUE))
    message("Directory '", output_dir, "' created successfully.")
  }

  records_by_sample <- split(ebi_antibiogram_table, ebi_antibiogram_table$biosample_id)

  for (biosample in names(records_by_sample)) {
    output_records <- list()

    for (entry in 1:nrow(records_by_sample[[biosample]])) {
      output_records[[entry]] <- list(
        "antibioticName" = list(
          value = records_by_sample[[biosample]]$antibiotic_name[entry],
          iri = NULL
        ),
        "astStandard" = list(
          value = records_by_sample[[biosample]]$ast_standard[entry],
          iri = NULL
        ),
        "breakpointVersion" = list(
          value = breakpoint_version,
          iri = NULL
        ),
        "laboratoryTypingMethod" = list(
          value = records_by_sample[[biosample]]$laboratory_typing_method[entry],
          iri = NULL
        ),
        "measurement" = list(
          value = records_by_sample[[biosample]]$measurement[entry],
          iri = "null"
        ),
        "measurementUnits" = list(
          value = records_by_sample[[biosample]]$measurement_units[entry],
          iri = NULL
        ),
        "measurementSign" = list(
          value = records_by_sample[[biosample]]$measurement_sign[entry],
          iri = NULL
        ),
        "resistancePhenotype" = list(
          value = records_by_sample[[biosample]]$resistance_phenotype[entry],
          iri = NULL
        ),
        "platform" = list(
          value = records_by_sample[[biosample]]$platform[entry],
          iri = NULL
        )
      )
    }

    biosample_amr_record <- list(
      accession = biosample,
      data = list(
        list(
          domain = domain,
          webinSubmissionAccountId = submission_account,
          type = "AMR",
          schema = NULL,
          content = output_records
        )
      )
    )

    json_outfile <- file.path(output_dir, paste0(biosample, ".json"))

    safe_execute(write_json(
      biosample_amr_record,
      json_outfile,
      pretty = TRUE,
      auto_unbox = TRUE,
      null = "null"
    ))
  }
}

# Helper functions
safe_execute <- function(expr) {
  tryCatch(
    {
      expr
    },
    error = function(e) {
      message("Error in executing command: ", e$message)
      return(NULL)
    }
  )
}
