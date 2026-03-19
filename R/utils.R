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

AMRgen_env <- new.env()

utils::globalVariables(c(
  ".",
  "..count..",
  ".data",
  ".estimate",
  ".estimator",
  ".metric",
  ".row_id",
  "2.5 %",
  "97.5 %",
  "ab",
  "ab_code",
  "ab_col",
  "ab_name",
  "AMR::clinical_breakpoints",
  "AMR_associated_publications",
  "amr_method",
  "amrfp_drugs",
  "amrfp_drugs_table",
  "antibiotic",
  "Antibiotic",
  "antibiotic_name",
  "antibiotics",
  "ast_standard",
  "Best_Hit_ARO",
  "binary_comb",
  "BioProject",
  "bioproject_acc",
  "biosample_acc",
  "biosample_id",
  "BioSample_ID",
  "breakpoint_R",
  "breakpoint_S",
  "CARD Short Name",
  "category",
  "ci.lower",
  "ci.lower.est",
  "ci.lower.ppv",
  "ci.upper",
  "ci.upper.est",
  "ci.upper.ppv",
  "ci_lower",
  "ci_upper",
  "class_col",
  "coalesce",
  "Collection Date",
  "collection_date",
  "colours_ppv",
  "combination_id",
  "count",
  "count_label",
  "Cut_Off",
  "Disk diffusio",
  "Disk diffusion (mm)",
  "disk",
  "disk_diffusion",
  "disk_dose",
  "disk_potency",
  "disk_raw",
  "drug_agent",
  "drug_agent_code",
  "drug_agent_name",
  "drug_class",
  "drug_class_agent",
  "drug_class_from_agent",
  "drug_class_internal",
  "drug_internal",
  "drug_name_raw",
  "drug_to_parse",
  "ecoff",
  "ecoff_disk",
  "ecoff_mic",
  "ecoff_MIC",
  "Element subtype",
  "Element symbol",
  "Element type",
  "element_subtype_col",
  "element_symbol",
  "element_symbol_col",
  "element_type_col",
  "est",
  "Estimate",
  "estimate",
  "eucast",
  "fill_col",
  "Freq",
  "Gene symbol",
  "gene",
  "gene_symbol_col",
  "geno_prediction",
  "get(assay)",
  "get(marker_col)",
  "group",
  "guideline",
  "Hierarchy node",
  "hierarchy_node",
  "I.denom",
  "I.n",
  "I.ppv",
  "I.se",
  "id",
  "import_amrfp_ebi",
  "index",
  "interp_raw",
  "is_screening",
  "Kleborate_Class",
  "kleborate_classes",
  "Lab ID",
  "Laboratory typing method",
  "Laboratory typing platform",
  "laboratory_typing_method",
  "laboratory_typing_platform",
  "marker",
  "marker.label",
  "marker_count",
  "marker_list",
  "Measuremen",
  "Measurement sign",
  "measurement",
  "Measurement",
  "measurement_sign",
  "measurement_units",
  "median",
  "Method",
  "method",
  "method_code",
  "metric",
  "MIC (mg/L)",
  "mic",
  "mic_raw",
  "Microorganism",
  "microorganism",
  "microorganism_code",
  "mics",
  "mo",
  "mutation",
  "na.omit",
  "Name",
  "name",
  "node",
  "NWT",
  "NWT.denom",
  "NWT.n",
  "NWT.ppv",
  "NWT.se",
  "NWT_pred",
  "Organism Name",
  "Organism",
  "organism",
  "outcome",
  "p",
  "pd",
  "perc",
  "pheno",
  "pheno_clsi_disk",
  "pheno_clsi_mic",
  "pheno_disk",
  "pheno_eucast_disk",
  "pheno_eucast_mic",
  "pheno_mic",
  "pheno_MIC",
  "pheno_provided",
  "pheno_screening",
  "pheno_trm",
  "phenotype",
  "phenotype-AMR_associated_publications",
  "phenotype-antibiotic_name",
  "phenotype-ast_standard",
  "phenotype-gen_measurement",
  "phenotype-laboratory_typing_method",
  "phenotype-organism",
  "phenotype-platform",
  "phenotype-resistance_phenotype",
  "platform",
  "point_size",
  "ppv",
  "Pr(>|z|)",
  "predNWT",
  "predR",
  "pval",
  "py",
  "py_run_string",
  "R",
  "R.denom",
  "R.n",
  "R.ppv",
  "R.se",
  "R_pred",
  "r_to_py",
  "reagent",
  "Resistance Mechanism",
  "Resistance phenotype",
  "resistance_phenotype",
  "RGI_DrugClassAgent",
  "rgi_drugs_table",
  "rgi_short_name_table",
  "row_idx",
  "row_pct",
  "Sample",
  "sample_name",
  "Scientifi",
  "Scientific name",
  "scientific_name",
  "se",
  "setNames",
  "sig_binary",
  "sir",
  "sir_exp",
  "sir_exp",
  "sir_inst",
  "sir_inst",
  "sir_interp",
  "sir_raw",
  "sir_value",
  "sirscan_codes",
  "site",
  "SNPs_in_Best_Hit_ARO",
  "solo",
  "Source",
  "Specimen date",
  "specimen_type",
  "spp_pheno",
  "standard",
  "strain",
  "subclass",
  "Subclass",
  "subclass_col",
  "subclass_to_parse",
  "subtype",
  "symbol",
  "test_method",
  "Testin",
  "Testing Date",
  "Testing standard",
  "testing_standard",
  "truth_value",
  "type",
  "type",
  "Type",
  "u",
  "user",
  "value",
  "variation type",
  "x"
))


# dummy function to stop note 'All declared Imports should be used'
ignore_unused_imports <- function() {
  rlang::sym
}

#' @importFrom readr read_tsv read_delim
process_input <- function(input) {
  if (is.character(input) && file.exists(input)) {
    tsv_ext <- paste(rep("tsv", 5), c("", ".gz", ".bz2", ".xz", ".zip"), sep = "", collapse = "|")
    txt_ext <- paste(rep("txt", 5), c("", ".gz", ".bz2", ".xz", ".zip"), sep = "", collapse = "|")
    csv_ext <- paste(rep("csv", 5), c("", ".gz", ".bz2", ".xz", ".zip"), sep = "", collapse = "|")
    if (grepl(tsv_ext, input)) {
      data <- readr::read_tsv(input)
    } else if (grepl(txt_ext, input)) {
      data <- readr::read_delim(input, delim = NULL)
    } else if (grepl(csv_ext, input)) {
      data <- readr::read_csv(input)
    }
  } else if (is.data.frame(input)) {
    # Check if the input is already a dataframe
    data <- input
  } else {
    # If the input is neither a file nor a dataframe, stop with an error
    stop("Input must be either a valid file path or a dataframe.")
  }
  # strip any leading hash (e.g. NCBI AST)
  data <- data %>% dplyr::rename_with(~ stringr::str_remove(.x, "#"))
  # Return the dataframe
  return(data)
}

#' @importFrom cli ansi_has_hyperlink_support
font_url <- function(url, txt = url) {
  if (tryCatch(isTRUE(ansi_has_hyperlink_support()), error = function(e) FALSE)) {
    paste0("\033]8;;", url, "\a", txt, "\033]8;;\a")
  } else {
    url
  }
}

#' @importFrom cli ansi_has_hyperlink_support
font_italic <- function(..., collapse = " ") {
  txt <- paste0(c(...), collapse = collapse)
  if (tryCatch(isTRUE(ansi_has_hyperlink_support()), error = function(e) FALSE)) {
    if (is.null(collapse)) {
      paste0("\033[3m", txt, "\033[23m", collapse = NULL)
    } else {
      paste0("\033[3m", txt, "\033[23m", collapse = "")
    }
  } else {
    txt
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

# TODO REMOVE THIS CODE FOR CRAN SUBMISSION
send_to_github <- function() {
  cli::cli_alert_info("Styling code using {.fn styler::style_pkg}...")
  st <- utils::capture.output(styler::style_pkg(style = styler::tidyverse_style))

  cli::cli_alert_info("Documenting code using {.fn devtools::document}...")
  doc <- devtools::document(quiet = TRUE)

  cli::cli_alert_info("Checking code using {.fn devtools::check}...")
  ch <- devtools::check(quiet = TRUE)

  if (length(ch$errors) > 0 || length(ch$warnings) > 0) {
    print(ch)
    cli::cli_alert_danger("Errors, warnings, and notes must be fixed before pushing to GitHub. You're almost there!")
    return(invisible())
  }

  cli::cli_alert_success("All tests passed!")
  commit_msg <- readline("Your commit message: ")
  q <- utils::askYesNo("Ready to push to GitHub?", prompts = c("Yes", "No", "Cancel"))
  if (isTRUE(q)) {
    system2("git", args = "add .")
    system2("git", args = paste("commit -m '", commit_msg, "'"))
    system2("git", args = "push")
    cli::cli_alert_success("Pushed to GitHub.")
  } else {
    cli::cli_alert_danger("Cancelled.")
  }
}
