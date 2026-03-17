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

#' Import and process antimicrobial phenotype data from a generic format
#'
#' This function attempts to import antibiotic susceptibility testing (AST) data in long-form antibiogram format (one row per sample and test), suitable for downstream use with AMRgen analysis functions. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be compressed) and parses relevant columns (antibiotic names, species names, MIC or disk data, S/I/R calls) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF).
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in long-form antibiogram format (one row per sample and test). This might be a file containing the content of an EBI/NCBI AST dataset previously processed using [import_ebi_ast()], [import_ncbi_ast()], or [import_ast()] functions, or files with a similar format/structure but with different column names.
#' @param sample_col (optional, default `"id"`) String indicating the name of the input data column that provides the sample name. If the `rename_cols` parameter is set to `TRUE`, this column will be renamed as `id`.
#' @param species_col (optional, default `"species"`) String indicating the name of the input data column that provides a species name. If provided, this column will be converted to micro-organism class `mo` via [AMR::as.mo()]. If the `rename_cols` parameter is set to `TRUE`, this column will also be renamed as `spp_pheno`. If interpretation is switched on, this column will be used to identify the appropriate  breakpoints for interpretation of each row in the data table.
#' @param species (optional, default `NULL`) String indicating the name of the single species to which all samples belong. Use this if you want to interpret assay measurements but the input file does not contain a column indicating the species for each sample (called `species` or a name specified by the `species_col` parameter).
#' @param ab_col (optional, default `"drug_agent"`) String indicating the name of the input data column that provides a drug name. If provided, this column will be converted to antibiotic class `ab` via [AMR::as.ab()]. If the `rename_cols` parameter is set to `TRUE`, this column will also be renamed as `drug_agent`. If interpretation is switched on, this column will be used to identify the appropriate breakpoints for interpretation of each row in the data table.
#' @param ab (optional, default `NULL`) String indicating the name of a single antibiotic to use for phenotype interpretation. Use this if you want to interpret assay measurements but the input file does not contain a column indicating the drug for each sample (called `drug_agent` or a name specified by the `ab_col` parameter).
#' @param mic_col (optional, default `"mic"`) String indicating the name of the input data column that provides MIC measurements. If provided, this column will be converted to MIC class `mic` via [AMR::as.mic()]. If the `rename_cols` parameter is set to `TRUE`, this column will also be renamed as `mic`. If interpretation is switched on, the MIC values will be interpreted against clinical breakpoints.
#' @param disk_col (optional, default `"disk"`) String indicating the name of the input data column that provides disk diffusion zone measurements. If provided, this column will be converted to disk diffusion class `disk` via [AMR::as.disk()]. If the `rename_cols` parameter is set to `TRUE`, this column will also be renamed as `disk`. If interpretation is switched on, the zone values will be interpreted against clinical breakpoints.
#' @param pheno_cols (optional, default `c("ecoff", "pheno_eucast", "pheno_clsi", "pheno_provided")`) String, or vector of strings, indicating the name/s of the input data column/s that provide interpreted phenotypes (with column values encoded as `S/I/R` or `WT/NWT`). If provided, these columns will be converted to class `sir` via [AMR::as.sir()].
#' @param method_col (optional, default `"method"`) String indicating the name of the input data column that indicates the testing method used (e.g. `"broth dilution"`, `"disk diffusion"`). If the `rename_cols` parameter is set to `TRUE`, this column will also be renamed as `method`.
#' @param platform_col (optional, default `"platform"`) String indicating the name of the input data column that indicates the testing platform used (e.g. `"Vitek"`, `"Sensititre"`). If the `rename_cols` parameter is set to `TRUE`, this column will also be renamed as `platform`.
#' @param source_col (optional, default `"source"`) String indicating the name of the input data column that indicates the source of the dataset (e.g. `"BioProject"`, `"PMID"`). If the `rename_cols` parameter is set to `TRUE`, this column will also be renamed as `source`.
#' @param guideline_col (optional, default `"guideline"`) String indicating the name of the input data column that indicates the guideline used for testing (e.g. `"EUCAST"`, `"CLSI"`). If the `rename_cols` parameter is set to `TRUE`, this column will also be renamed as `guideline`.
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @param rename_cols A logical value (default is `TRUE`). If `TRUE`, the function will rename the provided columns (specified by `sample_col`, `ab_col`, `mic_col`, `disk_col`, `species_col`) to the default names expected by AMRgen functions (`id`, `drug_agent`, `mic`, `disk`, `spp_pheno`), to match those output by the other [import_ast()] functions.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of mutate relocate
#' @importFrom rlang is_string :=
#' @return A data frame with the processed AST data, including additional columns:
#' @export
#' @examples
#' \dontrun{
#' # import and process AST data from EBI, write formatted data to file for later use
#' pheno <- import_ebi_ast("EBI_AMR_data.csv.gz")
#' write_tsv(pheno,
#'   file = "EBI_AMR_data_processed.tsv.gz",
#'   interpret_eucast = TRUE, interpret_ecoff = TRUE
#' )
#'
#' # read stored data and format the columns to the correct classes
#' pheno <- format_ast("EBI_AMR_data_processed.tsv.gz")
#'
#' # read in unprocessed E. coli AST data from non-standard format and interpret
#' pheno <- format_ast("AMR_data.tsv",
#'   sample_col = "STRAIN", species = "E. coli",
#'   ab_col = "Antibiotic", mic_col = "MIC (mg/L)",
#'   interpret_eucast = TRUE, interpret_ecoff = TRUE
#' )
#' }
format_ast <- function(input,
                       sample_col = "id",
                       species = NULL, # convert with as.mo
                       species_col = "spp_pheno", # convert with as.mo
                       ab = NULL, # convert with as.ab
                       ab_col = "drug_agent", # convert with as.ab
                       mic_col = "mic", # convert with as.mic
                       disk_col = "disk", # convert with as.disk
                       pheno_cols = c("ecoff", "pheno_eucast", "pheno_clsi", "pheno_provided"), # convert with as.sir
                       method_col = "method",
                       platform_col = "platform",
                       source_col = "source",
                       guideline_col = "guideline",
                       interpret_eucast = FALSE,
                       interpret_clsi = FALSE,
                       interpret_ecoff = FALSE,
                       rename_cols = TRUE) {
  ast <- process_input(input)

  if (!is.null(species)) {
    cat(paste("Adding new micro-organism column 'spp_pheno' (class 'mo') with constant value", species, "\n"))
    ast$spp_pheno <- as.mo(species)
  }

  if (!is.null(species_col)) {
    if (species_col %in% colnames(ast)) {
      ast <- ast %>% mutate(!!sym(species_col) := as.mo(!!sym(species_col)))
      cat(paste("Parsing column", species_col, "as micro-organism (class 'mo')\n"))
      if (rename_cols & species_col != "spp_pheno") {
        ast <- ast %>% rename(spp_pheno = !!sym(species_col))
        cat(paste("Renaming column", species_col, "to standard name 'spp_pheno'\n"))
      }
      if ((interpret_ecoff | interpret_eucast | interpret_clsi) & !("spp_pheno" %in% colnames(ast))) {
        ast <- ast %>% mutate(spp_pheno = !!sym(species_col)) # we need a column named spp_pheno for interpretation
      }
    } else {
      cat(paste("Could not find species_col", species_col, "in input table\n"))
    }
  }

  if (!is.null(ab)) {
    cat(paste("Adding new antibiotic column 'drug_agent' (class 'ab') with constant value", ab, "\n"))
    ast$drug_agent <- as.ab(ab)
  }

  if (!is.null(ab_col)) {
    if (ab_col %in% colnames(ast)) {
      ast <- ast %>% mutate(!!sym(ab_col) := as.ab(!!sym(ab_col)))
      cat(paste("Parsing column", ab_col, "as antibiotic (class 'ab')\n"))
      if (rename_cols & ab_col != "drug_agent") {
        ast <- ast %>% rename(drug_agent = !!sym(ab_col))
        cat(paste("Renaming column", ab_col, "to standard name 'drug_agent'\n"))
      }
    } else {
      cat(paste("Could not find ab_col", ab_col, "in input table\n"))
    }
  }

  if (!is.null(mic_col)) {
    if (mic_col %in% colnames(ast)) {
      ast <- ast %>% mutate(!!sym(mic_col) := as.mic(!!sym(mic_col)))
      cat(paste("Parsing column", mic_col, "as class 'mic'\n"))
      if (rename_cols & mic_col != "mic") {
        ast <- ast %>% rename(mic = !!sym(mic_col))
        cat(paste("Renaming column", mic_col, "to standard name 'mic'\n"))
      }
    } else {
      cat(paste("Could not find mic_col", mic_col, "in input table\n"))
    }
  }

  if (!is.null(disk_col)) {
    if (disk_col %in% colnames(ast)) {
      ast <- ast %>% mutate(!!sym(disk_col) := as.disk(!!sym(disk_col)))
      cat(paste("Parsing column", disk_col, "as class 'disk'\n"))
      if (rename_cols & disk_col != "disk") {
        ast <- ast %>% rename(disk = !!sym(disk_col))
        cat(paste("Renaming column", disk_col, "to standard name 'disk'\n"))
      }
    } else {
      cat(paste("Could not find disk_col", disk_col, "in input table\n"))
    }
  }

  if (!is.null(pheno_cols)) {
    if (is_string(pheno_cols)) { # single column
      if (pheno_cols %in% colnames(ast)) {
        ast <- ast %>% mutate(!!sym(pheno_cols) := as.sir(!!sym(pheno_cols)))
        cat(paste("Parsing column", pheno_cols, "as class 'sir'\n"))
      } else {
        cat(paste("Could not find pheno_cols", pheno_cols, "in input table\n"))
      }
    } else {
      for (pheno_col in pheno_cols) {
        if (!is.null(pheno_col)) {
          if (pheno_col %in% colnames(ast)) {
            ast <- ast %>% mutate(!!sym(pheno_col) := as.sir(!!sym(pheno_col)))
            cat(paste("Parsing column", pheno_col, "as class 'sir'\n"))
          } else {
            cat(paste("Could not find pheno_col", pheno_col, "in input table\n"))
          }
        }
      }
    }
  }

  if (!is.null(method_col)) {
    if (method_col %in% colnames(ast)) {
      if (rename_cols & method_col != "method") {
        ast <- ast %>% rename(method = !!sym(method_col))
        cat(paste("Renaming column", method_col, "to standard name 'method'\n"))
      }
    } else {
      cat(paste("Could not find method_col", method_col, "in input table\n"))
    }
  }

  if (!is.null(platform_col)) {
    if (platform_col %in% colnames(ast)) {
      if (rename_cols & platform_col != "platform") {
        ast <- ast %>% rename(platform = !!sym(platform_col))
        cat(paste("Renaming column", platform_col, "to standard name 'platform'\n"))
      }
    } else {
      cat(paste("Could not find platform_col", platform_col, "in input table\n"))
    }
  }

  if (!is.null(guideline_col)) {
    if (guideline_col %in% colnames(ast)) {
      if (rename_cols & guideline_col != "guideline") {
        ast <- ast %>% rename(guideline = !!sym(guideline_col))
        cat(paste("Renaming column", guideline_col, "to standard name 'guideline'\n"))
      }
    } else {
      cat(paste("Could not find guideline_col", guideline_col, "in input table\n"))
    }
  }

  if (!is.null(source_col)) {
    if (source_col %in% colnames(ast)) {
      if (rename_cols & source_col != "source") {
        ast <- ast %>% rename(source = !!sym(source_col))
        cat(paste("Renaming column", source_col, "to standard name 'source'\n"))
      }
    } else {
      cat(paste("Could not find source_col", source_col, "in input table\n"))
    }
  }

  ast <- interpret_ast(ast,
    interpret_ecoff = interpret_ecoff,
    interpret_eucast = interpret_eucast,
    interpret_clsi = interpret_clsi,
    species = species, # note column spp_pheno will be used otherwise
    ab = ab
  )

  if (rename_cols) {
    ast <- ast %>% rename(id = !!sym(sample_col))
  }

  ast <- ast %>% relocate(any_of(c(
    "id", "drug_agent", "mic", "disk", ab_col, mic_col, disk_col, pheno_cols,
    "pheno_eucast", "pheno_clsi", "ecoff", "guideline",
    "method", "platform", "source", "pheno_provided", "spp_pheno", species_col
  )))

  return(ast)
}


#' Interpret antimicrobial susceptibility phenotype data in a standard format tibble
#'
#' This function applies human EUCAST or CLSI breakpoints, and/or ECOFF, to interpret antimicrobial susceptibility testing (AST) data.
#' @param ast A tibble containing the phenotype measures in standard AMRgen format, as output by [import_ast()]. It must contain assay measurements in columns `mic` (class `mic`) and/or `disk` (class `disk`). Interpretation requires an organism (column `spp_pheno` of class `mo`, or a single value passed via the `species` parameter) and an antibiotic (column `drug_agent` of class `ab`, or a single value passed via the `ab` parameter).
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the `spp_pheno` field in the input file will be assumed to specify the species for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the `drug_agent` field in the input file will be assumed to specify the antibiotic for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr across coalesce mutate
#' @return A copy of the input table, with additional columns:
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as WT/NWT), based on the MIC or disk diffusion data.
#' - `spp_pheno`: The species phenotype, formatted using the [as.mo()] function (either taken from field `spp_pheno` in the input table, or the single value specified by `species` parameter).
#' - `drug_agent`: The antibiotic used in the test, formatted using the [as.ab()] function (either taken from field `drug_agent` in the input table, or the single value specified by `ab` parameter).
#' @export
#' @examples
#' # import without re-interpreting resistance
#' pheno <- import_ncbi_ast(ecoli_ast_raw)
#' head(pheno)
#'
#' # interpret phenotypes
#' pheno <- interpret_ast(pheno)
#'
#' \dontrun{
#' pheno <- read_csv("AST.csv") %>%
#'   # convert antibiotic field to 'drug_agent' of class 'ab'
#'   mutate(drug_agent = as.ab(antibiotic)) %>%
#'   mutate(mic = paste0(sign, MIC)) %>%
#'   mutate(mic = as.mic(mic)) # create a single 'mic' column of class 'mic'
#'
#' pheno <- interpret_ast(pheno, species = "Escherichia coli")
#' }
interpret_ast <- function(ast, interpret_ecoff = TRUE, interpret_eucast = TRUE, interpret_clsi = TRUE, species = NULL, ab = NULL) {
  if (interpret_ecoff | interpret_eucast | interpret_clsi) {
    # check we have species
    if (!is.null(species)) {
      cat(paste0("Interpreting all data as species: ", mo_name(as.mo(species)), "\n"))
      if ("spp_pheno" %in% colnames(ast)) {
        if (length(unique(ast$spp_pheno)) > 1) {
          cat("Warning: ignoring 'spp_pheno' column in input table, which contains multiple species:\n")
          cat(paste(unique(ast$spp_pheno), collapse = ", "))
          cat("\n")
        } else if (as.mo(unique(ast$spp_pheno)) != as.mo(species)) {
          cat(paste0("Warning: ignoring 'spp_pheno' column in input table, which indicates a different species: ", unique(ast$spp_pheno), "\n"))
        }
      }
      ast <- ast %>% mutate(spp_pheno = as.mo(species))
    } else if (!("spp_pheno" %in% colnames(ast))) {
      stop("Warning: Could not interpret data, need to provide 'species' parameter or 'spp_pheno' column in input table")
    }

    # check we have antibiotic
    if (!is.null(ab)) {
      cat(paste0("Interpreting all data as drug: ", ab_name(as.ab(ab)), "\n"))
      if ("drug_agent" %in% colnames(ast)) {
        if (length(unique(ast$drug_agent)) > 1) {
          cat("Warning: ignoring 'drug_agent' column in input table, which contains multiple drugs:\n")
          cat(paste(unique(ast$drug_agent), collapse = ", "))
          cat("\n")
        } else if (as.ab(unique(ast$drug_agent)) != as.ab(ab)) {
          cat(paste0("Warning: ignoring 'drug_agent' column in input table, which indicates a different drug: ", unique(ast$drug_agent), "\n"))
        }
      }
      ast <- ast %>% mutate(drug_agent = as.ab(ab))
    } else if (!("drug_agent" %in% colnames(ast))) {
      stop("Warning: Could not interpret data, need to provide 'species' parameter or 'drug_agent' column in input table")
    }

    # interpret data
    if (interpret_eucast) {
      if ("mic" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(
            across(
              where(is.mic),
              ~ as.sir(.x,
                mo = "spp_pheno", ab = "drug_agent",
                guideline = "EUCAST", capped_mic_handling = "conservative"
              ),
              .names = "pheno_eucast_mic"
            )
          )
      }
      if ("disk" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(
            across(
              where(is.disk),
              ~ as.sir(.x,
                mo = "spp_pheno", ab = "drug_agent",
                guideline = "EUCAST"
              ),
              .names = "pheno_eucast_disk"
            )
          )
      }
      if (("pheno_eucast_mic" %in% colnames(ast)) & ("pheno_eucast_disk" %in% colnames(ast))) {
        ast <- ast %>%
          mutate(pheno_eucast = coalesce(pheno_eucast_mic, pheno_eucast_disk))
      } else if ("pheno_eucast_mic" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_eucast = pheno_eucast_mic)
      } else if ("pheno_eucast_disk" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_eucast = pheno_eucast_disk)
      }
    }
    if (interpret_clsi) {
      if ("mic" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(
            across(
              where(is.mic),
              ~ as.sir(.x,
                mo = "spp_pheno", ab = "drug_agent",
                guideline = "CLSI", capped_mic_handling = "conservative"
              ),
              .names = "pheno_clsi_mic"
            )
          )
      }
      if ("disk" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(
            across(
              where(is.disk),
              ~ as.sir(.x,
                mo = "spp_pheno", ab = "drug_agent",
                guideline = "CLSI"
              ),
              .names = "pheno_clsi_disk"
            )
          )
      }
      if (("pheno_clsi_mic" %in% colnames(ast)) & ("pheno_clsi_disk" %in% colnames(ast))) {
        ast <- ast %>%
          mutate(pheno_clsi = coalesce(pheno_clsi_mic, pheno_clsi_disk))
      } else if ("pheno_clsi_mic" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_clsi = pheno_clsi_mic)
      } else if ("pheno_eucast_disk" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_clsi_disk = pheno_clsi_disk)
      }
    }
    if (interpret_ecoff) {
      if ("mic" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(
            across(
              where(is.mic),
              ~ as.sir(.x,
                mo = "spp_pheno", ab = "drug_agent",
                guideline = "EUCAST", breakpoint_type = "ECOFF",
                capped_mic_handling = "conservative"
              ),
              .names = "ecoff_mic"
            )
          )
      }
      if ("disk" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(
            across(
              where(is.disk),
              ~ as.sir(.x,
                mo = "spp_pheno", ab = "drug_agent",
                guideline = "EUCAST", breakpoint_type = "ECOFF"
              ),
              .names = "ecoff_disk"
            )
          )
      }
      if (("ecoff_mic" %in% colnames(ast)) & ("ecoff_disk" %in% colnames(ast))) {
        ast <- ast %>%
          mutate(ecoff = coalesce(ecoff_mic, ecoff_disk))
      } else if ("ecoff_mic" %in% colnames(ast)) {
        ast <- ast %>% rename(ecoff = ecoff_mic)
      } else if ("ecoff_disk" %in% colnames(ast)) {
        ast <- ast %>% rename(ecoff = ecoff_disk)
      }
    }
  }
  return(ast)
}

#' Import and process antimicrobial susceptibility phenotype data from the NCBI AST browser
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset, processes the data, and optionally interprets the results based on MIC or disk diffusion data. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be compressed) and parses relevant columns (antibiotic names, species names, MIC or disk data) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF).
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in NCBI antibiogram format. These files can be downloaded from NCBI AST browser, e.g. <https://www.ncbi.nlm.nih.gov/pathogens/ast#Pseudomonas%20aeruginosa>
#' @param sample_col A string indicating the name of the column with sample identifiers. If `NULL`, assume this is `"BioSample"`.
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the field `Scientific name` will be assumed to specify the species for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the field `Antibiotic` will be assumed to specify the antibiotic for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @param source (optional) A single string to record as the source of these data points, e.g. "NCBI_browser". By default, the field `BioProject` will be used to indicate the source for each row in the input file, but if this is missing or you want to override it with a single value for all samples, you may provide a source name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of coalesce if_else mutate relocate rename
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The biological sample identifier (renamed from `BioSample` or specified column).
#' - `spp_pheno`: The species phenotype, formatted using the [as.mo()] function.
#' - `drug_agent`: The antibiotic used in the test, formatted using the [as.ab()] function.
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the [as.mic()] function.
#' - `disk`: The disk diffusion measurement (in mm), formatted using the [as.disk()] function.
#' - `method`: The AST method (e.g., `"broth dilution"`, `"disk diffusion"`, `"Etest"`, `"agar dilution"`). Expected values are based on the NCBI/EBI antibiogram specification. Note that NCBI allows `"MIC"` as a synonym for `"broth dilution"` but this function will convert `"MIC"` to `"broth dilution"` on import.
#' - `platform`: The AST platform/instrument (e.g., `"Vitek"`, `"Phoenix"`, `"Sensititre"`).
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay.
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as WT/NWT), based on the MIC or disk diffusion data.
#' - `pheno_provided`: The original phenotype interpretation provided in the input file.
#' - `source`: The source of each data point (renamed from `BioProject` in the input file, or replaced with a single value passed in as the `source` parameter).
#' @export
#' @examples
#' # small example E. coli AST data from NCBI
#' ecoli_ast_raw
#'
#' # import without re-interpreting resistance
#' pheno <- import_ncbi_ast(ecoli_ast_raw)
#' head(pheno)
#'
#' # import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
#' pheno <- import_ncbi_ast(ecoli_ast_raw, interpret_eucast = TRUE, interpret_ecoff = TRUE)
#' head(pheno)
import_ncbi_ast <- function(input, sample_col = "BioSample", source = NULL, species = NULL, ab = NULL,
                            interpret_eucast = FALSE, interpret_clsi = FALSE, interpret_ecoff = FALSE) {
  ast <- process_input(input)

  # find id column
  if (!is.null(sample_col)) {
    if (sample_col %in% colnames(ast)) {
      ast <- ast %>% rename(id = any_of(sample_col))
    } else {
      stop(paste("Invalid column name:", sample_col))
    }
  } else {
    stop("Please specify the column containing sample identifiers, via parameter 'sample_col'")
  }

  if ("Laboratory typing platform" %in% colnames(ast)) {
    ast <- ast %>% mutate(platform = `Laboratory typing platform`)
  } else {
    cat("Warning: Expected AST platform column 'Laboratory typing platform' not found in input\n")
  }

  # parse method column
  # Note: NCBI antibiogram spec lists "MIC" as a synonym for "broth dilution"
  if ("Laboratory typing method" %in% colnames(ast)) {
    ast <- ast %>%
      mutate(method = `Laboratory typing method`) %>%
      mutate(method = if_else(!is.na(method) & method == "MIC", "broth dilution", method))
  } else {
    cat("Warning: Expected AST method column 'Laboratory typing method' not found in input\n")
    # if no method field, guess from MIC/disk column values
    if ("MIC (mg/L)" %in% colnames(ast) & "Disk diffusion (mm)" %in% colnames(ast)) {
      ast <- ast %>% mutate(method = case_when(
        !is.na(`MIC (mg/L)`) ~ "broth dilution",
        !is.na(`Disk diffusion (mm)`) ~ "disk diffusion",
        TRUE ~ NA
      ))
    }
    # sometimes etest may be indicated in lab typing platform, or method
    if ("Laboratory typing platform" %in% colnames(ast)) {
      ast <- ast %>% mutate(method = case_when(
        is.na(`Laboratory typing platform`) ~ method,
        grepl("etest", tolower(`Laboratory typing platform`)) ~ "Etest",
        grepl("e-test", tolower(`Laboratory typing platform`)) ~ "Etest",
        TRUE ~ method
      ))
    }
    if ("Laboratory typing method version or reagent" %in% colnames(ast)) {
      ast <- ast %>% mutate(method = case_when(
        is.na(`Laboratory typing method version or reagent`) ~ method,
        grepl("etest", tolower(`Laboratory typing method version or reagent`)) ~ "Etest",
        grepl("e-test", tolower(`Laboratory typing method version or reagent`)) ~ "Etest",
        TRUE ~ method
      ))
    }
  }

  # parse guideline column
  if ("Testing standard" %in% colnames(ast)) {
    ast <- ast %>% mutate(guideline = `Testing standard`)
  } else {
    cat("Warning: Expected AST standard column 'Testing standard' not found in input\n")
  }

  # parse disk column
  if ("Disk diffusion (mm)" %in% colnames(ast)) {
    ast <- ast %>% mutate(disk = as.disk(`Disk diffusion (mm)`))
  } else {
    cat("Warning: Expected column 'Disk diffusion (mm)' not found in input\n")
  }

  # parse mic column
  if ("MIC (mg/L)" %in% colnames(ast)) {
    if ("Measurement sign" %in% colnames(ast)) {
      ast <- ast %>%
        mutate(mic = paste0(`Measurement sign`, `MIC (mg/L)`)) %>%
        mutate(mic = gsub("==", "", mic))
    } else {
      ast <- ast %>% mutate(mic = `MIC (mg/L)`)
      cat("Warning: Expected column 'Measurement sign' not found in input, be careful interpreting MIC\n")
    }
    ast <- ast %>% mutate(mic = as.mic(mic))
  } else {
    cat("Warning: Expected column 'MIC (mg/L)' not found in input\n")
  }

  # parse antibiotic column
  if ("Antibiotic" %in% colnames(ast)) {
    ast <- ast %>% mutate(drug_agent = as.ab(Antibiotic))
  } else {
    stop("Expected column 'Antibiotic' not found in input.")
  }

  # parse species column
  if ("Scientific name" %in% colnames(ast)) {
    ast <- ast %>% mutate(spp_pheno = as.mo(`Scientific name`))
  } else {
    cat("Warning: Expected species column 'Scientific name' not found in input\n")
  }

  # parse bioproject column as source
  if ("BioProject" %in% colnames(ast)) {
    ast <- ast %>% mutate(source = BioProject)
  } else {
    if (!is.null(source)) {
      ast <- ast %>% mutate(source = source)
      cat(paste0("Setting source to user-provided value: ", source, "\n"))
    } else {
      cat("Warning: Expected column 'BioProject' not found in input\n")
    }
  }

  # parse phenotype SIR column
  if ("Resistance phenotype" %in% colnames(ast)) {
    ast <- ast %>%
      mutate(pheno_provided = as.sir(`Resistance phenotype`))
    if ("intermediate" %in% ast$`Resistance phenotype`) {
      message('Manually setting `pheno_provided` to "I" where `Resistance phenotype` was "intermediate"')
      ast <- ast %>%
        mutate(pheno_provided = if_else(`Resistance phenotype` == "intermediate",
          "I",
          `Resistance phenotype`
        ))
    }
  } else {
    cat("Warning: Expected phenotype SIR column 'Resistance phenotype' not found in input\n")
  }

  ast <- interpret_ast(ast, interpret_ecoff = interpret_ecoff, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, species = species, ab = ab)

  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", "pheno_eucast", "pheno_clsi", "ecoff", "guideline", "method", "platform", "source", "pheno_provided", "spp_pheno")))

  return(ast)
}

#' Import and process antimicrobial phenotype data retrieved from NCBI BioSamples
#'
#' This function will import antibiotic susceptibility testing (AST) data suitable for downstream use with AMRgen analysis functions. The expected input is phenotype data retrieved from NCBI BioSample database via the function [download_ncbi_ast()].
#' Note that files downloaded from the NCBI AST web browser <https://www.ncbi.nlm.nih.gov/pathogens/ast> are formatted differently and can be imported with the function [import_ncbi_ast()].
#' @param input A string representing the input dataframe, or a path to an input file, to be processed.
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of mutate relocate
#' @importFrom rlang is_string :=
#' @return A data frame with the processed AST data, including additional columns:
#' @export
#' @examples
#' \dontrun{
#' # Download Klebsiella quasipneumoniae data, filter to amikacin
#' ast <- download_ncbi_ast(
#'   "Klebsiella quasipneumoniae",
#'   antibiotic = "amikacin"
#' )
#'
#' # Reformat to simplify use with AMRgen functions
#' ast <- import_ncbi_biosample(ast, interpret_eucast = T)
#' }
import_ncbi_biosample <- function(input,
                                  interpret_eucast = FALSE,
                                  interpret_clsi = FALSE,
                                  interpret_ecoff = FALSE) {
  ast <- process_input(input)

  # Note: NCBI antibiogram spec lists "MIC" as a synonym for "broth dilution"
  # EUtils/API format has a single Measurement column; method indicates MIC vs disk diffusion
  ast <- ast %>%
    mutate(mic = if_else(`Laboratory typing method` == "MIC",
      paste0(`Measurement sign`, Measurement),
      NA
    )) %>%
    mutate(disk = if_else(`Laboratory typing method` == "disk diffusion",
      paste0(`Measurement sign`, Measurement),
      NA
    )) %>%
    mutate(`Laboratory typing method` = if_else(
      !is.na(`Laboratory typing method`) & `Laboratory typing method` == "MIC",
      "broth dilution",
      `Laboratory typing method`
    )) %>%
    mutate(pheno_provided = if_else(`Resistance phenotype` == "intermediate",
      "I",
      `Resistance phenotype`
    )) %>%
    format_ast(
      sample_col = "id",
      species_col = "organism",
      ab_col = "Antibiotic",
      pheno_cols = "pheno_provided",
      method_col = "Laboratory typing method",
      platform_col = "Laboratory typing platform",
      source_col = "BioProject",
      guideline_col = "Testing standard",
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff
    )

  return(ast)
}

#' Import and process antimicrobial susceptibility phenotype data from the EBI AMR web portal
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset that has been downloaded from the EBI AMR portal website (<https://www.ebi.ac.uk/amr/data/?view=experiments>)
#' Data downloaded from the EBI AMR Portal FTP site (<ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/>), either directly or via the function [download_ebi()], is formatted differently and can instead be processed using the [import_ebi_ast_ftp()] function.
#'
#' This function will process the data, and optionally interpret the results based on MIC or disk diffusion data. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be compressed) and parses relevant columns (antibiotic names, species names, MIC or disk measures) into suitable classes using the AMR package. It optionally can use the AMR package to interpret MIC/disk zone measures into susceptibility phenotypes (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF).
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in EBI antibiogram format. These files can be downloaded from the EBI AMR browser, e.g. <https://www.ebi.ac.uk/amr/data/?view=experiments>.
#' @param sample_col A string indicating the name of the column with sample identifiers. If `NULL`, assume this is `"phenotype-BioSample_ID"`.
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the field `phenotype-organism` will be assumed to specify the species for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the field `phenotype-antibiotic_name` will be assumed to specify the antibiotic for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @param source (optional) A single string to record as the source of these data points, e.g. "EBI_browser". By default, the field `phenotype-AMR_associated_publications` will be used to indicate the source for each row in the input file, but if this is missing or you want to override it with a single value for all samples, you may provide a source name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of coalesce if_else mutate relocate rename
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The biological sample identifier (renamed from `phenotype-BioSample_ID` or specified column).
#' - `spp_pheno`: The species phenotype, formatted using the [as.mo()] function.
#' - `drug_agent`: The antibiotic used in the test, formatted using the [as.ab()] function.
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the [as.mic()] function.
#' - `disk`: The disk diffusion measurement (in mm), formatted using the [as.disk()] function.
#' - `method`: The AST method (e.g., `"broth dilution"`, `"disk diffusion"`, `"Etest"`, `"agar dilution"`). Expected values are based on the NCBI/EBI antibiogram specification.
#' - `platform`: The AST platform/instrument (e.g., `"Vitek"`, `"Phoenix"`, `"Sensititre"`).
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay.
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as WT/NWT), based on the MIC or disk diffusion data.
#' - `pheno_provided`: The original phenotype interpretation provided in the input file.
#' - `source`: The source of each data point (renamed from the publications field in the input file, or replaced with a single value passed in as the `source` parameter).
#' @export
#' @examples
#' \dontrun{
#' # import without re-interpreting resistance
#' pheno <- import_ebi_ast("EBI_AMR_data.csv.gz")
#' head(pheno)
#'
#' # import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
#' pheno <- import_ebi_ast("EBI_AMR_data.csv.gz", interpret_eucast = TRUE, interpret_ecoff = TRUE)
#' head(pheno)
#' }
import_ebi_ast <- function(input, sample_col = "phenotype-BioSample_ID", source = NULL, species = NULL, ab = NULL,
                           interpret_eucast = FALSE, interpret_clsi = FALSE, interpret_ecoff = FALSE) {
  ast <- process_input(input)

  # find id column
  if (!is.null(sample_col)) {
    if (sample_col %in% colnames(ast)) {
      ast <- ast %>% rename(id = any_of(sample_col))
    } else {
      stop(paste("Invalid column name:", sample_col))
    }
  } else {
    stop("Please specify the column containing sample identifiers, via parameter 'sample_col'")
  }

  # parse disk data
  if ("phenotype-gen_measurement" %in% colnames(ast)) {
    ast <- ast %>%
      mutate(disk = if_else(grepl("mm", `phenotype-gen_measurement`),
        as.disk(`phenotype-gen_measurement`), NA
      )) %>%
      mutate(disk = as.disk(disk))
  } else {
    ast <- ast %>% mutate(disk = as.disk(NA))
    cat("No disk data (units 'mm') found in input\n")
  }

  # parse mic data
  if ("phenotype-gen_measurement" %in% colnames(ast)) {
    ast <- ast %>%
      mutate(mic = if_else(grepl("mg/L", `phenotype-gen_measurement`),
        as.mic(`phenotype-gen_measurement`), NA
      )) %>%
      mutate(mic = as.mic(mic))
  } else {
    ast <- ast %>% mutate(mic = as.mic(NA))
    cat("No MIC data (units 'mg/L') found in input\n")
  }

  # parse antibiotic column
  if ("phenotype-antibiotic_name" %in% colnames(ast)) {
    ast <- ast %>% mutate(drug_agent = as.ab(`phenotype-antibiotic_name`))
  } else {
    stop("Expected drug name column 'phenotype-antibiotic_name' not found in input.")
  }

  # parse species column
  if ("phenotype-organism" %in% colnames(ast)) {
    ast <- ast %>% mutate(spp_pheno = as.mo(`phenotype-organism`))
  } else {
    cat("Warning: Expected species column 'phenotype-organism' not found in input\n")
  }

  if ("phenotype-platform" %in% colnames(ast)) {
    ast <- ast %>% mutate(platform = `phenotype-platform`)
  } else {
    cat("Warning: Expected AST platform column 'phenotype-platform' not found in input\n")
  }

  if ("phenotype-laboratory_typing_method" %in% colnames(ast)) {
    ast <- ast %>% mutate(method = `phenotype-laboratory_typing_method`)
  } else {
    cat("Warning: Expected AST method column 'phenotype-laboratory_typing_method' not found in input\n")
  }

  if ("phenotype-ast_standard" %in% colnames(ast)) {
    ast <- ast %>% mutate(guideline = `phenotype-ast_standard`)
  } else {
    cat("Warning: Expected AST standard column 'phenotype-ast_standard' not found in input\n")
  }

  if ("phenotype-AMR_associated_publications" %in% colnames(ast)) {
    ast <- ast %>% mutate(source = `phenotype-AMR_associated_publications`)
  } else {
    if (!is.null(source)) {
      ast <- ast %>% mutate(source = source)
    }
    cat("Warning: Expected pubmed ID column 'phenotype-AMR_associated_publications' not found in input\n")
  }

  if ("phenotype-resistance_phenotype" %in% colnames(ast)) {
    ast <- ast %>% mutate(pheno_provided = as.sir(`phenotype-resistance_phenotype`))
  } else {
    cat("Warning: Expected pubmed ID column 'phenotype-resistance_phenotype' not found in input\n")
  }

  ast <- interpret_ast(ast, interpret_ecoff = interpret_ecoff, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, species = species, ab = ab)

  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", "pheno_eucast", "pheno_clsi", "ecoff", "guideline", "method", "platform", "source", "pheno_provided", "spp_pheno")))

  return(ast)
}


#' Import and process antimicrobial phenotype data files retrieved from the EBI AMR portal FTP site
#'
#' This function will import antibiotic susceptibility testing (AST) data suitable for downstream use with AMRgen analysis functions. The expected input is phenotype data retrieved from the [EBI AMR Portal FTP site](ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) either directly or via the function [download_ebi()].
#' Note that files downloaded from the [EBI AMR Portal web browser](https://www.ebi.ac.uk/amr/data/?view=experiments) are formatted differently and can be imported with the function [import_ebi_ast()].
#' @param input A string representing the input dataframe, or a path to an input file, to be processed.
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`..
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of mutate relocate
#' @importFrom rlang is_string :=
#' @return A data frame with the processed AST data, including additional columns:
#' @export
#' @examples
#' \dontrun{
#' # download Salmonella phenotype data from EBI
#' pheno_salmonella <- download_ebi(genus = "Salmonella")
#'
#' # reformat to simplify use with AMRgen functions
#' pheno_salmonella <- import_ebi_ast_ftp(pheno_salmonella)
#' }
import_ebi_ast_ftp <- function(input,
                               interpret_eucast = FALSE,
                               interpret_clsi = FALSE,
                               interpret_ecoff = FALSE) {
  ast <- process_input(input)

  ast <- ast %>%
    mutate(mic = if_else(measurement_units == "mg/L",
      paste0(measurement_sign, measurement),
      NA
    )) %>%
    mutate(disk = if_else(laboratory_typing_method == "disk diffusion",
      paste0(measurement_sign, measurement),
      NA
    )) %>%
    mutate(pheno_provided = if_else(resistance_phenotype == "intermediate",
      "I",
      resistance_phenotype
    )) %>%
    format_ast(
      sample_col = "BioSample_ID",
      species_col = "species",
      ab_col = "antibiotic_name",
      pheno_cols = "pheno_provided",
      method_col = "laboratory_typing_method",
      platform_col = "platform",
      source_col = "AMR_associated_publications",
      guideline_col = "ast_standard",
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff
    )

  return(ast)
}




#' Import and process antimicrobial phenotype data exported from Vitek instruments
#'
#' This function imports antimicrobial susceptibility testing (AST) data from Vitek instrument output files (wide CSV format)
#' and converts it to the standardised long-format used by AMRgen.
#'
#' @param input A dataframe or path to a CSV file containing Vitek AST output data
#' @param sample_col String indicating the name of the column containing sample identifiers. Default: `"Lab ID"`
#' @param source Optional string value to record in the `source` column for all data points (e.g., dataset name or study identifier)
#' @param species Optional string indicating a single species to use for phenotype interpretation (otherwise this is inferred per-sample from the input)
#' @param ab Optional string indicating a single antibiotic to use for phenotype interpretation (otherwise this is inferred per-sample from the input)
#' @param instrument_guideline Optional string indicating the guideline used by the instrument for SIR interpretation (e.g., "EUCAST 2025", "CLSI 2025"), used to record the `guideline` in the output file. Default: `NULL`
#' @param use_expertized Use expertized SIR (`TRUE`, default) or instrument SIR (`FALSE`)
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the wildtype vs nonwildtype status for each observation based on the MIC values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @param include_dates Logical indicating whether to include `collection_date` and `testing_date` in output. Default: `TRUE`
#' @importFrom AMR as.ab as.mic as.mo as.sir
#' @importFrom dplyr across any_of case_when coalesce mutate relocate rename
#' @importFrom tidyr pivot_longer matches
#' @importFrom stringr str_match
#' @importFrom rlang sym
#' @return Standardised AST data frame
#' @export
import_vitek_ast <- function(input,
                             sample_col = "Lab ID",
                             source = NULL,
                             species = NULL,
                             ab = NULL,
                             instrument_guideline = NULL,
                             use_expertized = TRUE,
                             interpret_eucast = FALSE,
                             interpret_clsi = FALSE,
                             interpret_ecoff = FALSE,
                             include_dates = TRUE) {
  ast <- process_input(input)

  # Validate sample column exists
  if (!(sample_col %in% colnames(ast))) {
    stop(paste("Invalid column name:", sample_col))
  }

  all_cols <- colnames(ast)

  # Drop 'Family -...' columns (resistance mechanism descriptions in lab VITEK files)
  family_cols <- all_cols[grepl("^Family -", all_cols)]
  if (length(family_cols) > 0) {
    cat(paste0("Dropping ", length(family_cols), " 'Family -...' resistance mechanism columns\n"))
    ast <- ast %>% dplyr::select(-any_of(family_cols))
    all_cols <- colnames(ast)
  }

  # Identify antibiotic MIC columns (CODE-Name pattern, excluding -Other-*)
  mic_cols <- all_cols[grepl("^[A-Z0-9]+-", all_cols) &
    !grepl("-Other-Instrument$|-Other-Expertized$", all_cols)]

  if (length(mic_cols) == 0) {
    stop("No antibiotic columns found in Vitek format")
  }

  # Extract antibiotic codes and names from column headers
  ab_codes <- stringr::str_match(mic_cols, "^([A-Z0-9]+)-")[, 2]
  ab_names <- stringr::str_match(mic_cols, "^[A-Z0-9]+-(.*)$")[, 2]

  # Metadata columns to preserve
  metadata_cols <- c(
    "Lab ID", "Isolate Number", "Patient Name", "Patient ID",
    "Patient Location", "Specimen Type", "Specimen Source",
    "Collection Date", "Testing Date", "Organism Name",
    "Organism Code", "Bio Number", "Percent Probability",
    "ID Confidence", "Selected BP Infection Site"
  )
  metadata_cols <- metadata_cols[metadata_cols %in% all_cols]

  # Rename columns for pivoting: CODE-Name -> CODE_mic, CODE-Other-Instrument -> CODE_sir_inst
  rename_map <- c()
  for (i in seq_along(ab_codes)) {
    code <- ab_codes[i]
    mic_col <- mic_cols[i]
    sir_inst_col <- paste0(code, "-Other-Instrument")
    sir_exp_col <- paste0(code, "-Other-Expertized")

    rename_map[paste0(code, "_mic")] <- mic_col
    if (sir_inst_col %in% all_cols) {
      rename_map[paste0(code, "_sir_inst")] <- sir_inst_col
    }
    if (sir_exp_col %in% all_cols) {
      rename_map[paste0(code, "_sir_exp")] <- sir_exp_col
    }
  }

  # Create lookup table for ab_code -> ab_name
  ab_lookup <- setNames(ab_names, ab_codes)

  # Apply renames
  for (new_name in names(rename_map)) {
    old_name <- rename_map[new_name]
    if (old_name %in% colnames(ast)) {
      colnames(ast)[colnames(ast) == old_name] <- new_name
    }
  }

  # Convert all MIC and SIR columns to character to avoid type conflicts in pivot
  mic_col_names <- paste0(ab_codes, "_mic")
  sir_inst_col_names <- paste0(ab_codes, "_sir_inst")
  sir_exp_col_names <- paste0(ab_codes, "_sir_exp")
  all_ab_cols <- c(mic_col_names, sir_inst_col_names, sir_exp_col_names)
  all_ab_cols <- all_ab_cols[all_ab_cols %in% colnames(ast)]

  ast <- ast %>%
    mutate(across(any_of(all_ab_cols), ~ as.character(.x)))

  # Pivot to long format
  ast_long <- ast %>%
    tidyr::pivot_longer(
      cols = tidyr::matches("^[A-Z0-9]+_(mic|sir_inst|sir_exp)$"),
      names_to = c("ab_code", ".value"),
      names_pattern = "^([A-Z0-9]+)_(mic|sir_inst|sir_exp)$"
    )

  # Add antibiotic full name and parse special values
  # (-) = not tested; TRM = too resistant to measure; POS/NEG = qualitative screening
  ast_long <- ast_long %>%
    mutate(ab_name = ab_lookup[ab_code]) %>%
    mutate(method = "broth dilution") %>%
    mutate(is_screening = mic %in% c("POS", "NEG")) %>%
    mutate(pheno_screening = case_when(
      mic == "POS" ~ "R",
      mic == "NEG" ~ "S",
      TRUE ~ NA_character_
    )) %>%
    mutate(pheno_trm = if_else(!is.na(mic) & mic == "TRM", "R", NA_character_)) %>%
    mutate(mic = case_when(
      is.na(mic) | mic == "" | mic == "-" | mic == "(-)" ~ NA_character_,
      mic %in% c("NEG", "POS", "TRM") ~ NA_character_,
      TRUE ~ mic
    )) %>%
    mutate(mic = as.mic(mic))

  # Select SIR source based on use_expertized parameter
  if (use_expertized && "sir_exp" %in% colnames(ast_long)) {
    ast_long <- ast_long %>%
      mutate(pheno_provided = case_when(
        sir_exp %in% c("S", "I", "R") ~ sir_exp,
        TRUE ~ NA_character_
      ))
  } else if ("sir_inst" %in% colnames(ast_long)) {
    ast_long <- ast_long %>%
      mutate(pheno_provided = case_when(
        sir_inst %in% c("S", "I", "R") ~ sir_inst,
        TRUE ~ NA_character_
      ))
  } else {
    ast_long <- ast_long %>% mutate(pheno_provided = NA_character_)
  }
  # Merge in screening and TRM results where no SIR was available
  ast_long <- ast_long %>%
    mutate(pheno_provided = coalesce(pheno_provided, pheno_screening, pheno_trm)) %>%
    mutate(method = if_else(is_screening, "screening", method)) %>%
    mutate(pheno_provided = as.sir(pheno_provided))

  # Parse organism and antibiotic
  if ("Organism Name" %in% colnames(ast_long)) {
    ast_long <- ast_long %>% mutate(spp_pheno = as.mo(`Organism Name`))
  }

  # Try full name first, then code for drug_agent (name is more reliable)
  ast_long <- ast_long %>%
    mutate(drug_agent_name = as.ab(ab_name)) %>%
    mutate(drug_agent_code = as.ab(ab_code)) %>%
    mutate(drug_agent = coalesce(drug_agent_name, drug_agent_code))

  # Rename sample column and add standard columns
  ast_long <- ast_long %>%
    rename(id = !!sym(sample_col)) %>%
    mutate(platform = "Vitek") %>%
    mutate(disk = as.disk(NA))

  # Add guideline if specified
  if (!is.null(instrument_guideline)) {
    ast_long <- ast_long %>% mutate(guideline = instrument_guideline)
  }

  # Add source only if specified
  if (!is.null(source)) {
    ast_long <- ast_long %>% mutate(source = source)
  }

  # Add dates if requested
  if (include_dates) {
    if ("Collection Date" %in% colnames(ast_long)) {
      ast_long <- ast_long %>% rename(collection_date = `Collection Date`)
    }
    if ("Testing Date" %in% colnames(ast_long)) {
      ast_long <- ast_long %>% rename(testing_date = `Testing Date`)
    }
  }

  # Interpret phenotypes
  ast_long <- interpret_ast(ast_long,
    interpret_ecoff = interpret_ecoff,
    interpret_eucast = interpret_eucast,
    interpret_clsi = interpret_clsi,
    species = species, ab = ab
  )

  # Reorder columns
  ast_long <- ast_long %>%
    relocate(any_of(c(
      "id", "drug_agent", "mic", "disk",
      "pheno_eucast", "pheno_clsi", "ecoff",
      "guideline", "method", "platform", "source",
      "pheno_provided", "spp_pheno"
    )))

  return(ast_long)
}


#' Import and process antimicrobial phenotype data exported from MicroScan instruments
#'
#' This function imports antimicrobial susceptibility testing (AST) data from MicroScan instrument output files (wide CSV format)
#' and converts it to the standardised long-format used by AMRgen. Supports English,
#' Spanish, French, German, and Portuguese column names (auto-detected from metadata columns).
#'
#' @param input A dataframe or path to a CSV/TSV file containing MicroScan AST output data
#' @param sample_col String indicating the name of the column containing sample identifiers. Default: `"Lab ID"`
#' @param source Optional string value to record in the `source` column for all data points (e.g., dataset name or study identifier)
#' @param species Optional string indicating a single species to use for phenotype interpretation (otherwise this is inferred per-sample from the input)
#' @param ab Optional string indicating a single antibiotic to use for phenotype interpretation (otherwise this is inferred per-sample from the input)
#' @param instrument_guideline Optional string indicating the guideline used by the instrument for SIR interpretation (e.g., "EUCAST 2025", "CLSI 2025"), used to record the `guideline` in the output file. Default: `NULL`
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the wildtype vs nonwildtype status for each observation based on the MIC values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr across any_of case_when coalesce filter mutate relocate rename select
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_match
#' @importFrom rlang sym
#' @return Standardised AST data frame
#' @export
import_microscan_ast <- function(input,
                                 sample_col = NULL,
                                 source = NULL,
                                 species = NULL,
                                 ab = NULL,
                                 instrument_guideline = NULL,
                                 interpret_eucast = FALSE,
                                 interpret_clsi = FALSE,
                                 interpret_ecoff = FALSE) {
  # MicroScan files may have .txt extension but be comma-separated
  # Try process_input first, fall back to read_csv if only 1 column detected
  ast <- process_input(input)
  if (ncol(ast) <= 1 && is.character(input) && file.exists(input)) {
    cat("Re-reading as CSV (file appears to be comma-separated despite .txt extension)\n")
    ast <- readr::read_csv(input, show_col_types = FALSE)
    ast <- ast %>% dplyr::rename_with(~ stringr::str_remove(.x, "#"))
  }

  # Strip leading/trailing single quotes from column names (some MicroScan exports)
  colnames(ast) <- gsub("^'|'$", "", colnames(ast))

  # Strip leading/trailing single quotes from all character columns
  ast <- ast %>%
    mutate(across(where(is.character), ~ gsub("^'|'$", "", .x)))

  # Detect language and translate metadata column names + drug suffixes
  # Supported: Spanish, French, German, Portuguese
  lang_markers <- list(
    spanish = list(
      detect = "Muestra",
      meta = c(
        "Muestra" = "Sample", "Origen" = "Source",
        "Aislamiento" = "Isolation", "Microorganismo" = "Microorganism"
      ),
      mic_suffix = " CIM",
      sir_suffixes = c(" Interpretaci\u00f3n", " Interpretacion")
    ),
    french = list(
      detect = c("\u00c9chantillon", "Echantillon"),
      meta = c(
        "\u00c9chantillon" = "Sample", "Echantillon" = "Sample",
        "Origine" = "Source", "Isolement" = "Isolation",
        "Micro-organisme" = "Microorganism", "Microorganisme" = "Microorganism"
      ),
      mic_suffix = " CMI",
      sir_suffixes = c(" Interpr\u00e9tation", " Interpretation")
    ),
    german = list(
      detect = "Probe",
      meta = c(
        "Probe" = "Sample", "Herkunft" = "Source", "Quelle" = "Source",
        "Isolierung" = "Isolation", "Mikroorganismus" = "Microorganism"
      ),
      mic_suffix = " MHK",
      sir_suffixes = c(" Interpretation")
    ),
    portuguese = list(
      detect = "Amostra",
      meta = c(
        "Amostra" = "Sample", "Origem" = "Source",
        "Isolamento" = "Isolation", "Microrganismo" = "Microorganism",
        "Microorganismo" = "Microorganism"
      ),
      mic_suffix = " CIM",
      sir_suffixes = c(" Interpreta\u00e7\u00e3o", " Interpretacao")
    )
  )

  detected_lang <- NULL
  for (lang_name in names(lang_markers)) {
    markers <- lang_markers[[lang_name]]$detect
    if (any(markers %in% colnames(ast))) {
      detected_lang <- lang_name
      break
    }
  }

  if (!is.null(detected_lang)) {
    lang <- lang_markers[[detected_lang]]
    cat(paste0("Detected ", detected_lang, "-language MicroScan export, translating column names\n"))

    # Translate metadata column names
    for (orig in names(lang$meta)) {
      if (orig %in% colnames(ast)) {
        colnames(ast)[colnames(ast) == orig] <- lang$meta[orig]
      }
    }

    # Translate MIC suffix
    colnames(ast) <- gsub(
      paste0(gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", lang$mic_suffix), "$"),
      " MIC", colnames(ast)
    )

    # Translate Interpretation suffixes
    for (suf in lang$sir_suffixes) {
      colnames(ast) <- gsub(
        paste0(gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", suf), "$"),
        " Interpretation", colnames(ast)
      )
    }
  }

  # Determine sample column (after translation, non-English names should be "Sample")
  if (is.null(sample_col)) {
    sample_candidates <- c(
      "Sample", "Muestra", "\u00c9chantillon", "Echantillon",
      "Probe", "Amostra"
    )
    found <- sample_candidates[sample_candidates %in% colnames(ast)]
    if (length(found) > 0) {
      sample_col <- found[1]
    } else {
      stop("Could not auto-detect sample column. Please specify via sample_col parameter.")
    }
  }
  if (!(sample_col %in% colnames(ast))) {
    stop(paste("Invalid column name:", sample_col))
  }

  # Identify drug columns by ' MIC' suffix
  all_cols <- colnames(ast)
  mic_cols <- all_cols[grepl(" MIC$", all_cols)]
  if (length(mic_cols) == 0) {
    stop("No antibiotic MIC columns found in MicroScan format")
  }

  # Extract drug names from column headers
  drug_names <- gsub(" MIC$", "", mic_cols)

  # Rename columns for pivoting: "DrugName MIC" -> "drugname___mic", "DrugName Interpretation" -> "drugname___sir"
  for (drug in drug_names) {
    mic_col_name <- paste0(drug, " MIC")
    sir_col_name <- paste0(drug, " Interpretation")
    safe_drug <- gsub("[^A-Za-z0-9]", "_", drug)
    if (mic_col_name %in% colnames(ast)) {
      colnames(ast)[colnames(ast) == mic_col_name] <- paste0(safe_drug, "___mic")
    }
    if (sir_col_name %in% colnames(ast)) {
      colnames(ast)[colnames(ast) == sir_col_name] <- paste0(safe_drug, "___sir")
    }
  }

  # Convert all drug columns to character
  drug_pattern_cols <- colnames(ast)[grepl("___mic$|___sir$", colnames(ast))]
  ast <- ast %>% mutate(across(any_of(drug_pattern_cols), as.character))

  # Pivot to long format
  ast_long <- ast %>%
    tidyr::pivot_longer(
      cols = tidyr::matches("^.+___(mic|sir)$"),
      names_to = c("drug_name_raw", ".value"),
      names_pattern = "^(.+)___(mic|sir)$"
    )

  # Restore drug name from safe version
  drug_lookup <- setNames(drug_names, gsub("[^A-Za-z0-9]", "_", drug_names))
  ast_long <- ast_long %>%
    mutate(drug_name_raw = drug_lookup[drug_name_raw])

  # Clean SIR values: strip asterisks (R* -> R, S* -> S, I* -> I)
  if ("sir" %in% colnames(ast_long)) {
    ast_long <- ast_long %>%
      mutate(sir_raw = sir) %>%
      mutate(sir = gsub("\\*$", "", sir)) %>%
      mutate(pheno_provided = case_when(
        sir %in% c("S", "I", "R") ~ sir,
        TRUE ~ NA_character_
      )) %>%
      mutate(pheno_provided = as.sir(pheno_provided))
  } else {
    ast_long <- ast_long %>% mutate(pheno_provided = as.sir(NA))
  }

  # Handle special MIC values
  ast_long <- ast_long %>%
    mutate(mic = case_when(
      is.na(mic) | mic == "" | mic == "(-)" ~ NA_character_,
      mic == "N/R" ~ NA_character_,
      mic == "ESBL" ~ NA_character_,
      TRUE ~ mic
    )) %>%
    mutate(mic = as.mic(mic))

  # Parse organism
  if ("Microorganism" %in% colnames(ast_long)) {
    ast_long <- ast_long %>% mutate(spp_pheno = as.mo(Microorganism))
  }

  # Parse drug_agent (drug names may be in Spanish; as.ab supports this)
  ast_long <- ast_long %>%
    mutate(drug_agent = as.ab(drug_name_raw))

  # Rename sample column and add standard columns
  ast_long <- ast_long %>%
    rename(id = !!sym(sample_col)) %>%
    mutate(method = "broth dilution") %>%
    mutate(platform = "MicroScan") %>%
    mutate(disk = as.disk(NA))

  # Add guideline if specified
  if (!is.null(instrument_guideline)) {
    ast_long <- ast_long %>% mutate(guideline = instrument_guideline)
  }

  # Add source only if specified
  if (!is.null(source)) {
    ast_long <- ast_long %>% mutate(source = source)
  }

  # Interpret phenotypes
  ast_long <- interpret_ast(ast_long,
    interpret_ecoff = interpret_ecoff,
    interpret_eucast = interpret_eucast,
    interpret_clsi = interpret_clsi,
    species = species, ab = ab
  )

  # Reorder columns
  ast_long <- ast_long %>%
    relocate(any_of(c(
      "id", "drug_agent", "mic", "disk",
      "pheno_eucast", "pheno_clsi", "ecoff",
      "guideline", "method", "platform", "source",
      "pheno_provided", "spp_pheno"
    )))

  return(ast_long)
}


#' Import and process antimicrobial phenotype data exported from Sensititre instruments
#'
#' This function imports antimicrobial susceptibility testing (AST) data from Sensititre instrument output files
#' (tab- or comma-separated, optionally UTF-16LE encoded, no header row) and converts it to the standardised
#' long-format used by AMRgen.
#'
#' @param input Path to a Sensititre output text file (tab- or comma-separated)
#' @param source Optional string value to record in the `source` column for all data points (e.g., dataset name or study identifier)
#' @param species Optional string indicating a single species to use for phenotype interpretation (otherwise this is inferred per-sample from the input)
#' @param ab Optional string indicating a single antibiotic to use for phenotype interpretation (otherwise this is inferred per-sample from the input)
#' @param instrument_guideline Optional string indicating the guideline used by the instrument for SIR interpretation (e.g., "EUCAST 2025", "CLSI 2025"), used to record the `guideline` in the output file. Default: `NULL`
#' @param id_col Integer. Column index (1-based) of the sample identifier. Default is 7, which corresponds
#'   to the sample accession column in standard Sensititre exports. Adjust if your file uses a different
#'   column for the sample ID (e.g. set to 2 for the plate/batch identifier column).
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the wildtype vs nonwildtype status for each observation based on the MIC values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of bind_rows case_when filter mutate relocate
#' @return Standardised AST data frame
#' @export
import_sensititre_ast <- function(input,
                                  source = NULL,
                                  species = NULL,
                                  ab = NULL,
                                  instrument_guideline = NULL,
                                  id_col = 7,
                                  interpret_eucast = FALSE,
                                  interpret_clsi = FALSE,
                                  interpret_ecoff = FALSE) {
  # Sensititre files may be UTF-16LE encoded; tab- or comma-separated; no header row
  # Cannot use process_input() - needs custom reading
  if (!is.character(input) || !file.exists(input)) {
    stop("import_sensititre_ast requires a file path as input")
  }

  # Detect encoding from BOM: UTF-16LE files start with 0xFF 0xFE
  bom <- readBin(input, "raw", n = 2)
  is_utf16le <- length(bom) >= 2 && bom[1] == as.raw(0xFF) && bom[2] == as.raw(0xFE)

  raw_lines <- if (is_utf16le) {
    tryCatch(
      {
        con <- file(input, encoding = "UTF-16LE")
        lines <- readLines(con, warn = FALSE)
        close(con)
        lines
      },
      error = function(e) readLines(input, warn = FALSE)
    )
  } else {
    readLines(input, warn = FALSE)
  }

  # Remove BOM character (if still present) and empty lines
  raw_lines <- sub("^\ufeff", "", raw_lines)
  raw_lines <- raw_lines[nchar(trimws(raw_lines)) > 0]

  if (length(raw_lines) == 0) {
    stop("No data found in Sensititre file")
  }

  # Auto-detect separator: count tabs vs commas in the first data line
  first_line <- raw_lines[1]
  n_tabs <- nchar(first_line) - nchar(gsub("\t", "", first_line, fixed = TRUE))
  n_commas <- nchar(first_line) - nchar(gsub(",", "", first_line, fixed = TRUE))
  sep <- if (n_tabs >= n_commas) "\t" else ","

  cat(paste0("Reading ", length(raw_lines), " rows from Sensititre file\n"))

  # Parse each line, extract drug triplets
  all_rows <- list()
  for (row_idx in seq_along(raw_lines)) {
    fields <- strsplit(raw_lines[row_idx], sep, fixed = TRUE)[[1]]
    # Strip surrounding double-quotes from all fields
    fields <- gsub('^"|"$', "", fields)

    # Find the datetime field to locate where metadata ends.
    # Supports both ISO format (YYYY-MM-DD HH:MM:SS) and European/US format
    # (DD/MM/YYYY HH:MM or MM/DD/YYYY HH:MM) as produced by some Sensititre exports.
    ts_idx <- which(grepl(
      "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}|\\d{2}/\\d{2}/\\d{4} \\d{2}:\\d{2}",
      fields
    ))[1]
    if (is.na(ts_idx)) {
      cat(paste0("Warning: No timestamp found in row ", row_idx, ", skipping\n"))
      next
    }

    # Extract metadata
    sample_id <- if (length(fields) >= id_col) trimws(fields[id_col]) else NA_character_
    panel_code <- if (length(fields) >= 8) trimws(fields[8]) else NA_character_
    organism_code <- if (length(fields) >= 10) trimws(fields[10]) else NA_character_
    specimen <- if (length(fields) >= 12) trimws(fields[12]) else NA_character_
    timestamp <- trimws(fields[ts_idx])

    # Parse drug data starting after the timestamp
    # Fields are in triplets (drug_name, mic_value, interpretation) but
    # empty positions may have irregular numbers of blank fields.
    # Strategy: scan for non-blank fields that look like drug names (all uppercase letters),
    # then take the next two fields as MIC and interpretation.
    drug_start <- ts_idx + 1
    if (drug_start > length(fields)) next
    remaining <- fields[drug_start:length(fields)]

    i <- 1
    while (i <= length(remaining) - 2) {
      field_val <- trimws(remaining[i])
      # Skip blank/empty fields
      if (field_val == "" || field_val == " " || is.na(field_val)) {
        i <- i + 1
        next
      }
      # This should be a drug name (uppercase letters/digits, 3-6 chars)
      if (grepl("^[A-Z][A-Z0-9]{1,5}$", field_val)) {
        drug_abbrev <- field_val
        mic_raw <- trimws(remaining[i + 1])
        interp_raw <- trimws(remaining[i + 2])

        all_rows[[length(all_rows) + 1]] <- data.frame(
          id = sample_id,
          panel = panel_code,
          organism_code = organism_code,
          specimen = specimen,
          testing_date = timestamp,
          drug_name_raw = drug_abbrev,
          mic_raw = mic_raw,
          interp_raw = interp_raw,
          stringsAsFactors = FALSE
        )
        i <- i + 3
      } else {
        i <- i + 1
      }
    }
  }

  if (length(all_rows) == 0) {
    stop("No drug data found in Sensititre file")
  }


  ast_long <- dplyr::bind_rows(all_rows) %>% as_tibble()

  # Clean MIC values: remove leading/trailing spaces, normalise signs
  ast_long <- ast_long %>%
    mutate(mic_raw = trimws(mic_raw)) %>%
    mutate(mic_raw = gsub("\\s+", " ", mic_raw)) %>%
    mutate(mic = case_when(
      is.na(mic_raw) | mic_raw == "" | mic_raw == " " ~ NA_character_,
      # Convert " = 4" -> "4", "<= 8" -> "<=8", " > 16" -> ">16"
      grepl("^[<>= ]+", mic_raw) ~ gsub("\\s+", "", mic_raw),
      TRUE ~ mic_raw
    )) %>%
    # Remove bare "=" prefix (just means exact value)
    mutate(mic = gsub("^=", "", mic)) %>%
    # Convert European decimal comma to period (e.g. "<=0,06" -> "<=0.06")
    mutate(mic = gsub(",", ".", mic, fixed = TRUE)) %>%
    mutate(mic = as.mic(mic))

  # Map interpretation codes to S/I/R
  ast_long <- ast_long %>%
    mutate(pheno_provided = case_when(
      interp_raw == "SUSC" ~ "S",
      interp_raw == "RESIST" ~ "R",
      interp_raw == "INTER" ~ "I", # Intermediate
      interp_raw == "SUSIE" ~ "I", # EUCAST "susceptible, increased exposure"
      interp_raw == "SC" ~ "I", # SDD (susceptible dose-dependent)
      interp_raw %in% c("NOINTP", "NOIMAT") ~ NA_character_, # no interpretive criteria
      interp_raw == "INEV" ~ NA_character_, # inherent/intrinsic result, no breakpoint applicable
      TRUE ~ NA_character_
    )) %>%
    mutate(pheno_provided = as.sir(pheno_provided))

  # Parse drug_agent from abbreviation
  ast_long <- ast_long %>%
    mutate(drug_agent = as.ab(drug_name_raw))

  # Add standard columns
  ast_long <- ast_long %>%
    mutate(method = "broth dilution") %>%
    mutate(platform = "Sensititre") %>%
    mutate(disk = as.disk(NA))

  # Add guideline if specified
  if (!is.null(instrument_guideline)) {
    ast_long <- ast_long %>% mutate(guideline = instrument_guideline)
  }

  # Add source only if specified
  if (!is.null(source)) {
    ast_long <- ast_long %>% mutate(source = source)
  }

  # Parse organism if present
  if ("organism_code" %in% colnames(ast_long)) {
    ast_long <- ast_long %>%
      mutate(spp_pheno = as.mo(organism_code))
  }

  # Interpret phenotypes
  ast_long <- interpret_ast(ast_long,
    interpret_ecoff = interpret_ecoff,
    interpret_eucast = interpret_eucast,
    interpret_clsi = interpret_clsi,
    species = species, ab = ab
  )

  # Reorder columns
  ast_long <- ast_long %>%
    relocate(any_of(c(
      "id", "drug_agent", "mic", "disk",
      "pheno_eucast", "pheno_clsi", "ecoff",
      "guideline", "method", "platform", "source",
      "pheno_provided", "spp_pheno"
    )))

  return(ast_long)
}


#' Import and process antimicrobial phenotype data from WHONET files
#'
#' This function imports antimicrobial susceptibility testing (AST) data from WHONET software output files (wide CSV format)
#' and converts it to the standardised long-format used by AMRgen.
#'
#' @param input A dataframe or path to a CSV file containing WHONET AST output data
#' @param sample_col Column name for sample identifiers. If `NULL` (default), the function
#'   auto-detects from known WHONET column names: `"Identification number"`, `"Identification"`,
#'   `"laboratory"`, `"patient_id"`. Supply a column name explicitly to override.
#' @param source Optional string value to record in the `source` column for all data points (e.g., dataset name or study identifier)
#' @param species Optional string indicating a single species to use for phenotype interpretation (otherwise this is inferred per-sample from the input)
#' @param ab Optional string indicating a single antibiotic to use for phenotype interpretation (otherwise this is inferred per-sample from the input)
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the wildtype vs nonwildtype status for each observation based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @param include_patient_info Include patient demographic columns in output
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr across any_of case_when coalesce mutate relocate rename
#' @importFrom tidyr pivot_longer matches
#' @importFrom stringr str_match
#' @importFrom rlang sym
#' @return Standardised AST data frame
#' @examples
#' # Built-in AMR package WHONET example dataset (standard uppercase format)
#' result <- import_whonet_ast(AMR::WHONET)
#' head(result)
#'
#' \dontrun{
#' # WHONET file with lowercase columns and _SIR suffix (e.g. amc_nd20, amc_nd20_SIR)
#' # import_whonet_ast("path/to/whonet_export.csv", sample_col = "patient_id")
#'
#' # Include patient demographics in output
#' import_whonet_ast("path/to/whonet_export.csv",
#'   sample_col = "patient_id",
#'   include_patient_info = TRUE
#' )
#'
#' # Interpret against EUCAST breakpoints for a specific species
#' result_interp <- import_whonet_ast(AMR::WHONET,
#'   interpret_eucast = TRUE,
#'   species = "Escherichia coli"
#' )
#' }
#' @export
import_whonet_ast <- function(input,
                              sample_col = NULL,
                              source = NULL,
                              species = NULL,
                              ab = NULL,
                              interpret_eucast = FALSE,
                              interpret_clsi = FALSE,
                              interpret_ecoff = FALSE,
                              include_patient_info = FALSE) {
  ast <- process_input(input)

  # Auto-detect or validate sample column
  known_sample_cols <- c("Identification number", "Identification", "laboratory", "patient_id")
  if (is.null(sample_col)) {
    sample_col <- known_sample_cols[known_sample_cols %in% colnames(ast)][1]
    if (is.na(sample_col)) {
      stop(paste(
        "Could not auto-detect sample column. Please specify sample_col. Looked for:",
        paste(known_sample_cols, collapse = ", ")
      ))
    }
  } else if (!(sample_col %in% colnames(ast))) {
    stop(paste("Invalid column name:", sample_col))
  }

  # Identify antibiotic columns (CODE_METHOD pattern, e.g., AMP_ND10, CIP_ED5, amc_nd20)
  # Also handles WHONET variant where SIR interpretation columns have a _SIR suffix
  # (e.g. amc_nd20 = raw zone/MIC value, amc_nd20_SIR = S/I/R interpretation)
  all_cols <- colnames(ast)

  # Detect _SIR suffix format
  sir_suffix_cols <- all_cols[grepl("^[A-Za-z]{2,4}_[A-Za-z0-9]+_SIR$", all_cols,
    ignore.case = TRUE
  )]
  has_sir_suffix <- length(sir_suffix_cols) > 0

  if (has_sir_suffix) {
    # Raw measurement columns are the same names with _SIR stripped
    ab_cols_base <- sub("_SIR$", "", sir_suffix_cols, ignore.case = TRUE)
    ab_cols <- ab_cols_base[ab_cols_base %in% all_cols]
    if (length(ab_cols) == 0) {
      # Only SIR columns present (no raw measurement columns) — treat normally
      ab_cols <- sir_suffix_cols
      has_sir_suffix <- FALSE
    }
  } else {
    ab_cols <- all_cols[grepl("^[A-Za-z]{2,4}_[A-Za-z]{2}", all_cols, ignore.case = TRUE)]
  }

  if (length(ab_cols) == 0) {
    stop("No antibiotic columns found in WHONET format")
  }

  # Metadata columns to preserve
  metadata_cols <- c(
    "Identification number", "Specimen number", "Organism",
    "Country", "Laboratory", "Specimen date", "Specimen type",
    "Isolate number", "Organism type", "Date of data entry"
  )
  if (include_patient_info) {
    metadata_cols <- c(
      metadata_cols, "Last name", "First name", "Sex",
      "Age", "Age category", "Date of admission"
    )
  }
  metadata_cols <- metadata_cols[metadata_cols %in% all_cols]

  # Convert antibiotic columns to character to avoid type conflicts
  ast <- ast %>%
    mutate(across(any_of(c(ab_cols, sir_suffix_cols)), as.character))

  # Pivot to long format; when _SIR suffix columns are present, also extract their
  # SIR interpretations and join them in so downstream parsing can prefer them
  if (has_sir_suffix) {
    ast$.row_id <- seq_len(nrow(ast))

    ast_long <- ast %>%
      tidyr::pivot_longer(
        cols = any_of(ab_cols),
        names_to = "ab_col",
        values_to = "sir_value"
      )

    # Pivot the _SIR columns, renaming them to their base names for the join
    sir_interp_long <- ast %>%
      select(.row_id, any_of(sir_suffix_cols)) %>%
      dplyr::rename_with(
        ~ sub("_SIR$", "", ., ignore.case = TRUE),
        any_of(sir_suffix_cols)
      ) %>%
      tidyr::pivot_longer(
        cols = any_of(ab_cols),
        names_to = "ab_col",
        values_to = "sir_interp"
      )

    ast_long <- ast_long %>%
      dplyr::left_join(sir_interp_long, by = c(".row_id", "ab_col")) %>%
      dplyr::select(-.row_id) %>%
      dplyr::select(-any_of(sir_suffix_cols))
  } else {
    ast_long <- ast %>%
      tidyr::pivot_longer(
        cols = any_of(ab_cols),
        names_to = "ab_col",
        values_to = "sir_value"
      )
    ast_long$sir_interp <- NA_character_
  }

  # Parse antibiotic directly from column name (as.ab handles WHONET format)
  ast_long <- ast_long %>%
    mutate(drug_agent = as.ab(ab_col))


  # Parse method, platform, and guideline from column name suffix
  #
  # Format: [GUIDELINE/PLATFORM][METHOD][POTENCY] where:

  # Standard codes (guideline-based):
  #   - Guideline: N = CLSI/NCCLS, E = EUCAST, D = DIN
  #   - Method: D = Disk diffusion, M = MIC, E = Etest
  #   - Potency: disk strength in ug (for disk diffusion only)
  #   e.g., ND10 = CLSI disk 10ug, EM = EUCAST MIC, EE = EUCAST Etest
  #
  # Platform codes (instrument-specific):
  #   - Platform: V = Vitek, P = Phoenix, M = Microscan, S = Sensititre, K = Trek
  #   - Method: M = MIC, D = Disk, K = Kirby-Bauer disk
  #   e.g., VM = Vitek MIC, PD = Phoenix Disk
  #
  # Reference: https://whonet.org/WebDocs/WHONET%202.Laboratory%20configuration.html
  ast_long <- ast_long %>%
    mutate(method_code = toupper(stringr::str_match(ab_col, "(?i)^[A-Za-z]{2,4}_(.*)$")[, 2])) %>%
    # Parse guideline (for standard codes N/E/D)
    mutate(guideline = case_when(
      grepl("^N", method_code) ~ "CLSI",
      grepl("^E", method_code) ~ "EUCAST",
      grepl("^D[DM]", method_code) ~ "DIN", # DD or DM for DIN (not platform codes like VD)
      TRUE ~ NA_character_
    )) %>%
    # Parse platform (for platform codes V/P/M/S/K)
    mutate(platform = case_when(
      grepl("^V", method_code) ~ "Vitek",
      grepl("^P", method_code) ~ "Phoenix",
      grepl("^M[MD]", method_code) ~ "Microscan", # MM or MD
      grepl("^S", method_code) ~ "Sensititre",
      grepl("^K", method_code) ~ "Trek",
      TRUE ~ NA_character_
    )) %>%
    # Parse method (D=disk, M=MIC/broth dilution, E=Etest, K=Kirby-Bauer)
    mutate(method = case_when(
      # Standard disk diffusion with potency
      grepl("^[NED]D[0-9]", method_code) ~ "disk diffusion",
      # Platform disk
      grepl("^[VPMSK]D", method_code) ~ "disk diffusion",
      grepl("^[VPMSK]K", method_code) ~ "disk diffusion",
      # Standard MIC (broth dilution)
      grepl("^[NED]M", method_code) ~ "broth dilution",
      # Platform MIC (broth dilution)
      grepl("^[VPMSK]M", method_code) ~ "broth dilution",
      # Etest
      grepl("^[NED]E", method_code) ~ "Etest",
      TRUE ~ NA_character_
    )) %>%
    # Extract disk potency for disk diffusion
    mutate(disk_potency = case_when(
      grepl("^[NED]D[0-9]", method_code) ~ gsub("^[NED]D", "", method_code),
      TRUE ~ NA_character_
    ))

  # Parse SIR values — prefer the explicit _SIR column when present, then fall back
  # to reading S/I/R directly from the raw value column (standard WHONET format)
  ast_long <- ast_long %>%
    mutate(pheno_provided = case_when(
      !is.na(sir_interp) & sir_interp %in% c("S", "I", "R") ~ sir_interp,
      sir_value %in% c("S", "I", "R") ~ sir_value,
      TRUE ~ NA_character_
    )) %>%
    mutate(pheno_provided = as.sir(pheno_provided)) %>%
    dplyr::select(-sir_interp)

  # Extract numeric measurements from sir_value where applicable.
  # WHONET column values may contain MIC strings (e.g. "<=0.5", ">32") for
  # broth-dilution/Etest columns, or zone sizes (e.g. "24") for disk columns;
  # as.mic() / as.disk() silently return NA for non-parseable values (e.g. "S").
  ast_long <- ast_long %>%
    mutate(
      mic = as.mic(case_when(
        method %in% c("broth dilution", "Etest") ~ sir_value,
        TRUE ~ NA_character_
      )),
      disk = as.disk(case_when(
        method == "disk diffusion" ~ sir_value,
        TRUE ~ NA_character_
      ))
    )

  # Parse organism
  if ("Organism" %in% colnames(ast_long)) {
    ast_long <- ast_long %>% mutate(spp_pheno = as.mo(Organism))
  }

  # Rename sample column
  ast_long <- ast_long %>%
    rename(id = !!sym(sample_col))

  # Add source only if specified
  if (!is.null(source)) {
    ast_long <- ast_long %>% mutate(source = source)
  }

  # Rename date columns if present
  if ("Specimen date" %in% colnames(ast_long)) {
    ast_long <- ast_long %>% rename(collection_date = `Specimen date`)
  }

  # Interpret phenotypes
  ast_long <- interpret_ast(ast_long,
    interpret_ecoff = interpret_ecoff,
    interpret_eucast = interpret_eucast,
    interpret_clsi = interpret_clsi,
    species = species, ab = ab
  )

  # Drop patient demographic columns when not requested (case-insensitive to handle
  # WHONET files that use lowercase column names)
  if (!include_patient_info) {
    patient_cols_standard <- c(
      "Last name", "First name", "Sex", "Age",
      "Age category", "Date of admission"
    )
    cols_to_drop <- colnames(ast_long)[
      tolower(colnames(ast_long)) %in% tolower(patient_cols_standard)
    ]
    ast_long <- ast_long %>% dplyr::select(-any_of(cols_to_drop))
  }

  # Reorder columns
  ast_long <- ast_long %>%
    relocate(any_of(c(
      "id", "drug_agent", "mic", "disk",
      "pheno_eucast", "pheno_clsi", "ecoff",
      "guideline", "method", "platform", "disk_potency", "source",
      "pheno_provided", "spp_pheno"
    )))

  return(ast_long)
}


#' Import and process antimicrobial phenotype data exported from BD Phoenix instruments
#'
#' This function imports antimicrobial susceptibility testing (AST) data from BD Phoenix instrument
#' export files and converts them to the standardised long-format used by AMRgen.
#'
#' Phoenix instruments can export data in different file formats. This function handles both
#' named-header exports (e.g. tab-delimited CLSI reports with columns such as
#' `"Antimicrobial"`, `"MIC or Concentration"`, `"Final (SIR)"`) and
#' positional/headerless exports (XLS files with no column headers). Column detection is
#' automatic based on common Phoenix header name patterns, but any column can be overridden
#' by supplying its name or 1-based index via `drug_col`, `mic_col`, `sir_col`,
#' `sample_col`, and `species_col`.
#'
#' For headerless files (detected when no standard column names are recognised), the function
#' assumes the standard Phoenix positional layout: column 1 = sample ID, column 2 = species,
#' column 3 = drug, column 4 = MIC, column 5 = instrument interpretation, column 6 =
#' expert/final interpretation. Drug names are standardised using [AMR::as.ab()], which handles
#' abbreviated, full, and language-specific names (e.g. `"Fosfomycin mit G6P"`). Species
#' names are standardised using [AMR::as.mo()].
#'
#' For tab-delimited files, trailing non-data sections commonly appended by Phoenix (e.g.
#' `"Resistance Markers"`, `"Expert Triggered Rules"`) are automatically discarded.
#'
#' @param input Path to a Phoenix export file (XLS, XLSX, TXT, or TSV), or a data frame.
#' @param sample_col Name or 1-based column index of the sample ID column. If `NULL`
#'   (default), auto-detected from common names `"Sample"`, `"ID"`, `"Isolate"`).
#'   For single-isolate files with no sample column, the filename is used as the sample ID.
#' @param drug_col Name or 1-based column index of the drug/antimicrobial column. If `NULL`
#'   (default), auto-detected from names such as `"Antimicrobial"`, `"Drug"`, etc.
#' @param mic_col Name or 1-based column index of the MIC column. If NULL` (default),
#'   auto-detected from names such as `"MIC"`, `"MIC or Concentration"`, etc.
#' @param sir_col Name or 1-based column index of the S/I/R interpretation column. If `NULL`
#'   (default), auto-detected based on `use_expertized`: prefers `"Final (SIR)"`,
#'   `"Expert (SIR)"` when `TRUE`; prefers `"Interp"` when `FALSE`.
#' @param species_col Name or 1-based column index of the species/organism column. If `NULL`
#'   (default), auto-detected from names such as `"Organism"`, `"Species"`, etc.
#'   Ignored if `species` is supplied.
#' @param species Optional fixed species string applied to all rows (parsed via [AMR::as.mo()]).
#'   Overrides `species_col`. Required for files without an organism column.
#' @param ab Optional antibiotic filter/override passed to [interpret_ast()].
#' @param use_expertized If `TRUE` (default), prefer the expert/final interpretation column
#'   over raw instrument interpretation when both are present and `sir_col` is not specified.
#' @param source Optional source label recorded in the `source` output column.
#' @param instrument_guideline Optional guideline string recorded in the `guideline` column
#'   (e.g. `"CLSI 2018"`, `"EUCAST 2024"`).
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class `sir`.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the susceptibility phenotype (SIR) for each observation based on the MIC values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class `sir`.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will re-interpret the wildtype vs nonwildtype status for each observation based on the MIC values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class `sir` and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of filter if_else mutate relocate
#' @importFrom readr read_tsv
#' @return Standardised long-format AST data frame
#' @export
#' @examples
#' \dontrun{
#' # Named-header export (e.g. CLSI report, auto-detected)
#' pheno <- import_phoenix_ast("phoenix_clsi_report.txt",
#'   species = "Escherichia coli", instrument_guideline = "CLSI 2018"
#' )
#'
#' # Positional/headerless XLS export (auto-detected)
#' pheno <- import_phoenix_ast("phoenix_export.xls")
#'
#' # Override column detection when names differ from defaults
#' pheno <- import_phoenix_ast("phoenix_export.xlsx",
#'   drug_col = "Antibiotic", mic_col = "Result",
#'   sir_col = "SIR_final", sample_col = "IsolateID",
#'   species_col = "Organism"
#' )
#' }
import_phoenix_ast <- function(input,
                               sample_col = NULL,
                               drug_col = NULL,
                               mic_col = NULL,
                               sir_col = NULL,
                               species_col = NULL,
                               species = NULL,
                               ab = NULL,
                               use_expertized = TRUE,
                               source = NULL,
                               instrument_guideline = NULL,
                               interpret_eucast = FALSE,
                               interpret_clsi = FALSE,
                               interpret_ecoff = FALSE) {
  is_file <- is.character(input) && file.exists(input)
  ext <- if (is_file) tolower(tools::file_ext(input)) else NULL

  # -------------------------------------------------------------------
  # Helper: resolve a column reference (name or 1-based index)
  # -------------------------------------------------------------------
  .resolve_col <- function(col_ref, df) {
    if (is.null(col_ref)) {
      return(NULL)
    }
    if (is.numeric(col_ref)) {
      if (col_ref < 1 || col_ref > ncol(df)) {
        stop(paste("Column index", col_ref, "is out of range (", ncol(df), "columns)."))
      }
      return(colnames(df)[col_ref])
    }
    if (!col_ref %in% colnames(df)) {
      stop(paste0(
        "Column '", col_ref, "' not found.\nAvailable columns: ",
        paste(colnames(df), collapse = ", ")
      ))
    }
    col_ref
  }

  # -------------------------------------------------------------------
  # Helper: auto-detect a column by regex patterns (case-insensitive)
  # -------------------------------------------------------------------
  .find_col <- function(df, patterns) {
    cols_lower <- tolower(colnames(df))
    for (p in patterns) {
      idx <- grep(p, cols_lower)[1]
      if (!is.na(idx)) {
        return(colnames(df)[idx])
      }
    }
    NULL
  }

  # -------------------------------------------------------------------
  # Helper: normalise MIC strings
  #   - Unicode <= sign (\u2264) -> <=
  #   - Comma decimal separator -> period (safety for locale variants)
  #   - Strip combination denominators: ">32/2" -> ">32"
  #   - Empty / dash -> NA
  # -------------------------------------------------------------------
  .clean_mic <- function(x) {
    x <- trimws(as.character(x))
    x[x == "" | x == "-"] <- NA_character_
    x <- gsub("\u2264", "<=", x)
    x <- gsub(",", ".", x)
    x <- gsub("/[0-9.,\u2264<=]+$", "", x)
    x
  }

  # -------------------------------------------------------------------
  # 1. Read raw data
  # -------------------------------------------------------------------
  if (is_file && ext %in% c("xls", "xlsx")) {
    # Peek at column names to determine whether the file has a header row
    peek <- tryCatch(
      readxl::read_excel(input, col_names = TRUE, .name_repair = "minimal", n_max = 0),
      error = function(e) NULL
    )
    known_pats <- c(
      "antimicrobial", "antibiotic", "drug", "^mic$",
      "interp", "interpretation", "^sample$", "organism", "species"
    )
    has_header <- !is.null(peek) &&
      any(sapply(known_pats, function(p) any(grepl(p, tolower(colnames(peek))))))

    ast <- readxl::read_excel(input, col_names = has_header, .name_repair = "minimal")
    if (!has_header) colnames(ast) <- paste0("V", seq_len(ncol(ast)))
  } else if (is_file && ext %in% c("txt", "tsv")) {
    lines <- readLines(input, warn = FALSE)
    # Locate the header row (first line starting with a recognised field name)
    header_idx <- which(grepl("^(Antimicrobial|Drug|Sample|Antibiotic)",
      lines,
      ignore.case = TRUE
    ))[1]
    if (is.na(header_idx)) header_idx <- 1
    # Strip trailing non-data sections appended by Phoenix
    stop_pat <- "^\\s*$|^Resistance Markers|^Expert Triggered Rules|^Interpretation Criteria"
    stop_cands <- which(grepl(stop_pat, lines) & seq_along(lines) > header_idx)
    data_end <- if (length(stop_cands) > 0) min(stop_cands) - 1 else length(lines)
    ast <- readr::read_tsv(
      I(paste(lines[header_idx:data_end], collapse = "\n")),
      show_col_types = FALSE, name_repair = "minimal"
    ) %>% dplyr::filter(!is.na(.[[1]]) & trimws(.[[1]]) != "")
  } else if (is.data.frame(input)) {
    ast <- input
  } else {
    ast <- process_input(input)
  }

  # -------------------------------------------------------------------
  # 2. Resolve column names
  # -------------------------------------------------------------------

  # Apply any user-supplied overrides first
  r_sample <- .resolve_col(sample_col, ast)
  r_drug <- .resolve_col(drug_col, ast)
  r_mic <- .resolve_col(mic_col, ast)
  r_sir <- .resolve_col(sir_col, ast)
  r_species <- .resolve_col(species_col, ast)

  # Auto-detect columns not yet resolved
  if (is.null(r_drug)) {
    r_drug <- .find_col(ast, c("antimicrobial", "antibiotic", "^drug$", "^agent$", "^ab$"))
  }
  if (is.null(r_mic)) {
    r_mic <- .find_col(ast, c("^mic$", "mic or concentration", "concentration"))
  }
  if (is.null(r_sir)) {
    if (use_expertized) {
      r_sir <- .find_col(ast, c("final.*sir", "expert.*sir", "^expertized$"))
      if (is.null(r_sir)) {
        r_sir <- .find_col(ast, c("^interp$", "^interpretation$", "^sir$"))
      }
    } else {
      r_sir <- .find_col(ast, c("^interp$", "^interpretation$"))
      if (is.null(r_sir)) {
        r_sir <- .find_col(ast, c("final.*sir", "expert.*sir", "^sir$"))
      }
    }
  }
  if (is.null(r_sample)) {
    r_sample <- .find_col(ast, c("^sample$", "^id$", "^isolate$", "sample_id", "patient_id"))
  }
  if (is.null(r_species) && is.null(species)) {
    r_species <- .find_col(ast, c("^organism$", "^species$", "^spp$", "^pathogen$"))
  }

  # Positional fallback for headerless files (columns named V1, V2, ...)
  if (is.null(r_drug) && all(grepl("^V\\d+$", colnames(ast)))) {
    if (ncol(ast) >= 4) {
      message(
        "No standard column headers detected.\n",
        "Using Phoenix positional defaults: col 1 = sample, col 2 = species, ",
        "col 3 = drug, col 4 = MIC, col 5 = instrument SIR, col 6 = expert SIR.\n",
        "Supply 'drug_col', 'mic_col', 'sir_col' etc. to override."
      )
      if (is.null(r_sample)) r_sample <- colnames(ast)[1]
      if (is.null(r_species)) r_species <- if (ncol(ast) >= 2) colnames(ast)[2] else NULL
      r_drug <- colnames(ast)[3]
      r_mic <- colnames(ast)[4]
      r_sir <- if (use_expertized && ncol(ast) >= 6) colnames(ast)[6] else colnames(ast)[5]
    }
  }

  # Validate that all required columns were resolved
  if (is.null(r_drug)) {
    stop(
      "Could not auto-detect the drug/antimicrobial column.\n",
      "Available columns: ", paste(colnames(ast), collapse = ", "),
      "\nUse 'drug_col' to specify."
    )
  }
  if (is.null(r_mic)) {
    stop(
      "Could not auto-detect the MIC column.\n",
      "Available columns: ", paste(colnames(ast), collapse = ", "),
      "\nUse 'mic_col' to specify."
    )
  }
  if (is.null(r_sir)) {
    stop(
      "Could not auto-detect the S/I/R interpretation column.\n",
      "Available columns: ", paste(colnames(ast), collapse = ", "),
      "\nUse 'sir_col' to specify."
    )
  }

  # Handle missing sample column
  if (is.null(r_sample)) {
    if (is_file) {
      derived_id <- gsub("\\.[^.]*$", "", basename(input))
      message(
        "No sample ID column found -- using filename '", derived_id,
        "' as sample ID. Use 'sample_col' to specify the sample ID column."
      )
      ast[["__id__"]] <- derived_id
      r_sample <- "__id__"
    } else {
      message("No sample ID column found. Use 'sample_col' to specify.")
    }
  }

  if (is.null(r_species) && is.null(species)) {
    message(
      "No species column found and 'species' not supplied.\n",
      "Breakpoint interpretation will not be possible without organism information.\n",
      "Use 'species_col' or 'species' to provide organism information."
    )
  }

  # -------------------------------------------------------------------
  # 3. Map to standardised columns and apply AMR package classes
  # -------------------------------------------------------------------

  # "X" = no breakpoint defined in Phoenix -> NA
  sir_vals <- as.character(ast[[r_sir]])
  sir_vals[sir_vals %in% c("X", "x", "-") | is.na(sir_vals)] <- NA_character_

  ast <- ast %>%
    mutate(
      id             = if (!is.null(r_sample)) as.character(.data[[r_sample]]) else NA_character_,
      drug_agent     = as.ab(.data[[r_drug]]),
      mic            = as.mic(.clean_mic(as.character(.data[[r_mic]]))),
      pheno_provided = as.sir(sir_vals)
    )

  if (!is.null(species)) {
    ast <- ast %>% mutate(spp_pheno = as.mo(species))
  } else if (!is.null(r_species)) {
    ast <- ast %>% mutate(spp_pheno = as.mo(.data[[r_species]]))
  }

  # -------------------------------------------------------------------
  # 4. Common processing
  # -------------------------------------------------------------------
  ast <- ast %>%
    mutate(
      platform = "Phoenix",
      method = "broth dilution",
      disk = as.disk(NA)
    )

  if (!is.null(instrument_guideline)) ast <- ast %>% mutate(guideline = instrument_guideline)
  if (!is.null(source)) ast <- ast %>% mutate(source = source)

  ast <- interpret_ast(ast,
    interpret_ecoff = interpret_ecoff,
    interpret_eucast = interpret_eucast,
    interpret_clsi = interpret_clsi,
    species = species, ab = ab
  )

  ast %>% relocate(any_of(c(
    "id", "drug_agent", "mic", "disk",
    "pheno_eucast", "pheno_clsi", "ecoff",
    "guideline", "method", "platform", "source",
    "pheno_provided", "spp_pheno"
  )))
}


#' Import and process antimicrobial phenotype data from common sources
#'
#' This function imports an antibiotic susceptibility testing datasets in formats exported by EBI, NCBI, WHOnet and several automated AST instruments (Vitek, Microscan, Sensititre, Phoenix). It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF).
#' @param input A string representing a dataframe, or a path to an input file, containing the phenotype data in a supported format. These files may be downloaded from public sources such as the [EBI AMR web browser](https://www.ebi.ac.uk/amr/data/?view=experiments), [EBI FTP site](ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/), or [NCBI browser](https://www.ncbi.nlm.nih.gov/pathogens/ast), or using the functions [download_ebi()], [download_ncbi_ast()], or [query_ncbi_bq_geno()]; or the files may be exported from supported AST instruments.
#' @param format A string indicating the format of the data: `"ebi"` (default), `"ebi_web"`, `"ebi_ftp"`, `"ncbi"`, `"ncbi_biosample"`, `"vitek"`, `"microscan"`, `"phoenix"`, `"sensititre"`, or `"whonet"`. This determines which importer function the data is passed on to for processing (see below).
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @param ... Format-specific arguments. See
#' - `"ebi"` : [import_ebi_ast()]
#' - `"ebi_web"` : [import_ebi_ast()]
#' - `"ebi_ftp"` :[import_ebi_ast_ftp()]
#' - `"ncbi"` : [import_ncbi_ast()]
#' - `"ncbi_biosample"` : [import_ncbi_biosample()]
#' - `"vitek"` : [import_vitek_ast()]
#' - `"microscan"` : [import_microscan_ast()]
#' - `"sensititre"` : [import_sensititre_ast()]
#' - `"phoenix"` : [import_phoenix_ast()]
#' - `"whonet"` : [import_whonet_ast()]
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The sample identifier (`character`).
#' - `spp_pheno`: The species phenotype, formatted using the [AMR::as.mo()] function (class `mo`).
#' - `drug_agent`: The antibiotic used in the test, formatted using the [AMR::as.ab()] function (class `ab`).
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the [AMR::as.mic()] function (class `mic`).
#' - `disk`: The disk diffusion measurement (in mm), formatted using the [AMR::as.disk()] function (class `disk`).
#' - `method`: The AST method (e.g., `"broth dilution"`, `"disk diffusion"`, `"Etest"`, `"agar dilution"`). Expected values are based on the NCBI/EBI antibiogram specification (`character`).
#' - `platform`: The AST platform/instrument (e.g., `"Vitek"`, `"Phoenix"`, `"Sensititre"`) (`character`).
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay (`character`).
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as `S/I/R`), based on the MIC or disk diffusion data (class `sir`).
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as `S/I/R`), based on the MIC or disk diffusion data (class `sir`).
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as `WT/NWT`), based on the MIC or disk diffusion data (class `sir`).
#' - `pheno_provided`: The original phenotype interpretation provided in the input file, formatted using [AMR::as.sir()] (class `sir`).
#' - `source`: The source of each data point (from the publications or bioproject field in the input file, or replaced with a single value passed in as the `source` parameter) (`character`).
#' @export
#' @examples
#' \dontrun{
#' # import NCBI data retrieved from Google Cloud, without re-interpreting resistance
#' head(staph_ast_ncbi_cloud_raw)
#' pheno <- import_pheno(staph_ast_ncbi_cloud_raw, format = "ncbi")
#'
#' # import NCBI data where biosample column has been renamed to 'id'
#' head(staph_ast_ncbi_raw)
#' import_pheno(staph_ast_ncbi_raw, "ncbi", sample_col = "id")
#'
#' # import NCBI data and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF)
#' head(ecoli_ast_raw)
#' pheno <- import_pheno(ecoli_ast_raw,
#'   format = "ncbi",
#'   interpret_eucast = TRUE, interpret_ecoff = TRUE
#' )
#'
#' # download Klebsiella quasipneumoniae phenotype data from NCBI BioSample
#' kquasi_raw_ncbi <- download_ncbi_ast("Klebsiella quasipneumoniae")
#' head(kquasi_raw_ncbi)
#' # import the data and interpret against EUCAST breakpoints
#' pheno <- import_pheno(kquasi_raw_ncbi,
#'   format = "ncbi_biosample",
#'   interpret_eucast = T
#' )
#'
#' # download Klebsiella quasipneumoniae phenotype data from EBI
#' kquasi_raw_ebi <- download_ebi(species = "Klebsiella quasipneumoniae")
#' head(kquasi_raw_ebi)
#' # import the data and interpret against ecoff
#' pheno <- import_pheno(kquasi_raw_ebi,
#'   format = "ebi_ftp",
#'   interpret_ecoff = TRUE
#' )
#'
#' # import Vitek data from file, with default parameters
#' pheno <- import_pheno("vitek_export.tsv",
#'   format = "vitek"
#' )
#'
#' # import Vitek data from file
#' # specify guideline that was used, remove dates, ignore expertized calls
#' pheno <- import_pheno("vitek_export.tsv",
#'   format = "vitek",
#'   instrument_guideline = "EUCAST 2025",
#'   use_expertized = FALSE,
#'   include_dates = FALSE
#' )
#' }
import_pheno <- function(input,
                         format = "ebi",
                         interpret_eucast = FALSE,
                         interpret_clsi = FALSE,
                         interpret_ecoff = FALSE,
                         ...) {
  fun <- get_pheno_importer(format)

  dots <- forward_args(fun, ...)

  rlang::exec(fun,
    input = input,
    interpret_eucast = interpret_eucast,
    interpret_clsi = interpret_clsi,
    interpret_ecoff = interpret_ecoff,
    !!!dots
  )
}

# registry of pheno data import functions, to be dispatched via import_pheno()
.pheno_importers <- list(
  ebi = "import_ebi_ast",
  ebi_web = "import_ebi_ast",
  ebi_ftp = "import_ebi_ast_ftp",
  ncbi = "import_ncbi_ast",
  ncbi_biosample = "import_ncbi_biosample",
  vitek = "import_vitek_ast",
  microscan = "import_microscan_ast",
  sensititre = "import_sensititre_ast",
  phoenix = "import_phoenix_ast",
  whonet = "import_whonet_ast"
)

# function to identify the right phenotype importer function to dispatch
# case-insensitive to user-supplied parameter string 'format'
get_pheno_importer <- function(format) {
  format <- tolower(format)

  fun_name <- .pheno_importers[[format]]

  if (is.null(fun_name)) {
    stop(
      "Unknown format: ", format,
      "\nAvailable formats: ",
      paste(names(.pheno_importers), collapse = ", "),
      call. = FALSE
    )
  }

  get(fun_name, envir = asNamespace(utils::packageName()))
}

# forward relevant arguments to the function selected for dispatch
# silently drop any user-provided arguments that are not relevant to the function
forward_args <- function(fun, ...) {
  dots <- list(...)
  allowed <- names(formals(fun))

  keep <- names(dots) %in% allowed

  if (any(!keep)) {
    warning(
      "Ignoring unsupported arguments: ",
      paste(names(dots)[!keep], collapse = ", "),
      call. = FALSE
    )
  }

  dots[keep]
}


#' Import and process antimicrobial phenotype data from common sources
#'
#' This function imports an antibiotic susceptibility testing (AST) datasets in formats exported by EBI, NCBI, WHOnet and several automated AST instruments (Vitek, Microscan, Sensititre). It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be compressed) and parses relevant columns (antibiotic names, species names, MIC or disk data) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF). If expected columns are not found warnings will be given, and interpretation may not be possible.
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data a supported format. These files may be downloaded from public sources such as the EBI AMR web browser (<https://www.ebi.ac.uk/amr/data/?view=experiments>), EBI FTP site (<ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/>), or NCBI browser (e.g. <https://www.ncbi.nlm.nih.gov/pathogens/ast#Pseudomonas%20aeruginosa>), or using the functions [download_ebi()] or [download_ncbi_ast()]; or the files may be exported from supported AST instruments.
#' @param format A string indicating the format of the data: `"ebi"` (default), `"ebi_web"`, `"ebi_ftp"`, `"ncbi"`, `"vitek"`, `"microscan"`, `"sensititre"`, or `"whonet"`. This determines whether the data is passed on to the [import_ebi_ast()] (ebi/ebi_web), [import_ebi_ast_ftp()] (ebi_ftp), [import_ncbi_ast()] (ncbi), [import_vitek_ast()] (vitek), [import_microscan_ast()] (microscan), [import_sensititre_ast()] (sensititre), or [import_whonet_ast()] (whonet) function to process.
#' @param interpret_eucast A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is `FALSE`). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is `FALSE`). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as `NWT` (nonwildtype) or `WT` (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the field `organism` in the input file will be assumed to specify the species for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the field `antibiotic` in the input file will be assumed to specify the antibiotic for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @param source (optional) A string to record as the source of these data points, e.g. `"EBI_browser"`. By default, the publications field (for EBI data) or BioProject field (for NCBI data) will be used to indicate the source for each row in the input file. For other input files, or to override this with a single value for all samples, you may provide a source name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of coalesce if_else mutate relocate rename
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The biosample identifier.
#' - `spp_pheno`: The species phenotype, formatted using the [as.mo()] function.
#' - `drug_agent`: The antibiotic used in the test, formatted using the [as.ab()] function.
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the [as.mic()] function.
#' - `disk`: The disk diffusion measurement (in mm), formatted using the [as.disk()] function.
#' - `method`: The AST method (e.g., `"broth dilution"`, `"disk diffusion"`, `"Etest"`, `"agar dilution"`). Expected values are based on the NCBI/EBI antibiogram specification.
#' - `platform`: The AST platform/instrument (e.g., `"Vitek"`, `"Phoenix"`, `"Sensititre"`).
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay.
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as `S/I/R`), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as `S/I/R`), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as `WT/NWT`), based on the MIC or disk diffusion data.
#' - `pheno_provided`: The original phenotype interpretation provided in the input file.
#' - `source`: The source of each data point (from the publications or bioproject field in the input file, or replaced with a single value passed in as the `source` parameter).
#' @export
#' @examples
#' # small example E. coli AST data from NCBI
#' ecoli_ast_raw
#'
#' # import without re-interpreting resistance
#' pheno <- import_ast(ecoli_ast_raw, format = "ncbi")
#' head(pheno)
#'
#' # import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
#' pheno <- import_ast(ecoli_ast_raw, format = "ncbi", interpret_eucast = TRUE, interpret_ecoff = TRUE)
#' head(pheno)
import_ast <- function(input, format = "ebi", interpret_eucast = FALSE,
                       interpret_clsi = FALSE, interpret_ecoff = FALSE,
                       species = NULL, ab = NULL, source = NULL) {
  if (format %in% c("ebi_web", "ebi")) {
    cat("Reading in as EBI AST format downloaded from the web portal\n")
    ast <- import_ebi_ast(input,
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff,
      species = species, ab = ab
    )
  }

  if (format == "ebi_ftp") {
    cat("Reading in as EBI AST format downloaded from the FTP portal\n")
    ast <- import_ebi_ast_ftp(input,
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff
    )
  }

  if (format == "ncbi") {
    cat("Reading in as NCBI AST format\n")
    ast <- import_ncbi_ast(input, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, interpret_ecoff = interpret_ecoff, species = species, ab = ab)
  }

  if (format == "vitek") {
    cat("Reading in as VITEK AST format\n")
    ast <- import_vitek_ast(input,
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff,
      species = species, ab = ab, source = source
    )
  }

  if (format == "microscan") {
    cat("Reading in as MicroScan AST format\n")
    ast <- import_microscan_ast(input,
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff,
      species = species, ab = ab, source = source
    )
  }

  if (format == "sensititre") {
    cat("Reading in as Sensititre AST format\n")
    ast <- import_sensititre_ast(input,
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff,
      species = species, ab = ab, source = source
    )
  }

  if (format == "whonet") {
    cat("Reading in as WHONET AST format\n")
    ast <- import_whonet_ast(input,
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff,
      species = species, ab = ab, source = source
    )
  }

  if (format == "phoenix") {
    cat("Reading in as BD Phoenix AST format\n")
    ast <- import_phoenix_ast(input,
      interpret_eucast = interpret_eucast,
      interpret_clsi = interpret_clsi,
      interpret_ecoff = interpret_ecoff,
      species = species, ab = ab, source = source
    )
  }

  if (!is.null(source)) {
    ast <- ast %>% mutate(source = source)
  }
  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", "pheno_eucast", "pheno_clsi", "ecoff", "guideline", "method", "platform", "source", "pheno_provided", "spp_pheno")))

  return(ast)
}


#' Import phenotype predictions from AMRrules output
#'
#' This imports antimicrobial phenotype predictions (S/I/R and WT/NWT) generated from genotype data using [AMRrules](https://www.amrrules.org).
#' @param input A string representing a dataframe, or a path to an input file, containing the AMRrules output "genome_summary" file, which should be a long-form TSV file with one row per sample and drug.
#' @param sample_col (optional, default `"sample"`) String indicating the name of the input data column that provides the sample name.
#' @param species_col (optional, default `"organism"`) String indicating the name of the input data column that provides a species name. If provided, this column will be converted to micro-organism class `mo` via [AMR::as.mo()]. If the `rename_cols` parameter is set to `TRUE`, this column will also be renamed as `spp_pheno`. If interpretation is switched on, this column will be used to identify the appropriate  breakpoints for interpretation of each row in the data table.
#' @param ab_col (optional, default `"drug"`) String indicating the name of the input data column that provides a drug name.
#' @param sir_col (optional, default `"clinical category"`) String indicating the name of the input data column that indicates the S/I/R prediction.
#' @param ecoff_col (optional, default `"phenotype"`) String indicating the name of the input data column that indicates the WT/NWT prediction.
#' @param method (optional, default `"genotyping"`) String indicating the value to record in a new `method` field added to the output table.
#' @param platform (optional, default `"AMRfinderplus + AMRrules"`) String indicating the value to record in a new `platform` field added to the output table.
#' @importFrom AMR as.ab as.mo as.sir
#' @importFrom dplyr any_of mutate relocate
#' @return A data frame with the processed AST data, including additional columns:
#' @export
#' @examples
#' \dontrun{
#' # import and process AST data from EBI, write formatted data to file for later use
#' predictions <- import_amrrules_predictions("Ecoli_genome_summary.tsv")
#' }
import_amrrules_predictions <- function(input,
                                        sample_col = "sample",
                                        species_col = "organism", # convert with as.mo
                                        ab_col = "drug", # convert with as.ab
                                        sir_col = "clinical category", # convert with as.sir
                                        ecoff_col = "phenotype", # convert with as.sir
                                        method = "genotyping",
                                        platform = "AMRfinderplus + AMRrules") {
  intable <- process_input(input)

  if (!is.null(sample_col)) {
    if (sample_col %in% colnames(intable)) {
      intable <- intable %>%
        rename(id = !!sym(sample_col))
    }
  }

  if (!is.null(species_col)) {
    if (species_col %in% colnames(intable)) {
      intable <- intable %>%
        mutate(!!sym(species_col) := gsub("s__", "", !!sym(species_col))) %>%
        mutate(spp_pheno = as.mo(!!sym(species_col)))
    }
  }
  if (!is.null(ab_col)) {
    if (ab_col %in% colnames(intable)) {
      intable <- intable %>%
        mutate(drug_agent = as.ab(!!sym(ab_col)))
    }
  }
  if (!is.null(method)) {
    intable <- intable %>% mutate(method = method)
  }
  if (!is.null(platform)) {
    intable <- intable %>% mutate(platform = platform)
  }

  if (!is.null(sir_col)) {
    if (sir_col %in% colnames(intable)) {
      intable <- intable %>%
        mutate(!!sym(sir_col) := if_else(is.na(!!sym(sir_col)), "S", !!sym(sir_col))) %>%
        mutate(!!sym(sir_col) := as.sir(!!sym(sir_col)))
      cat(paste("Parsing column", sir_col, "as class 'sir'\n"))
    } else {
      cat(paste("Could not find SIR prediction column", sir_col, "in input table\n"))
    }
  }
  if (!is.null(ecoff_col)) {
    if (ecoff_col %in% colnames(intable)) {
      intable <- intable %>%
        mutate(!!sym(ecoff_col) := case_when(
          !!sym(ecoff_col) == "nonwildtype" ~ "NWT",
          !!sym(ecoff_col) == "wildtype" ~ "WT",
          !!sym(ecoff_col) == "-" ~ "WT",
          TRUE ~ NA
        )) %>%
        mutate(!!sym(ecoff_col) := as.sir(!!sym(ecoff_col)))
      cat(paste("Parsing column", ecoff_col, "as class 'sir'\n"))
    } else {
      cat(paste("Could not find WT/NWT phenotype prediction column", ecoff_col, "in input table\n"))
    }
  }

  intable <- intable %>%
    relocate(any_of(c("id", "drug_agent", sir_col, ecoff_col, "mo")))

  return(intable)
}
