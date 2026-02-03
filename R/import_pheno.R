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

#' Import and Process AST Data from an NCBI File
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset, processes the data, and optionally interprets the results based on MIC or disk diffusion data. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be compressed) and parses relevant columns (antibiotic names, species names, MIC or disk data) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF). If expected columns are not found warnings will be given, and interpretation may not be possible.
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in NCBI antibiogram format. These files can be downloaded from NCBI AST browser, e.g. https://www.ncbi.nlm.nih.gov/pathogens/ast#Pseudomonas%20aeruginosa
#' @param sample_col A string indicating the name of the column with sample identifiers. If `NULL`, assume this is 'BioSample'.
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the field 'Scientific name' will be assumed to specify the species for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the field 'Antibiotic' will be assumed to specify the antibiotic for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @param source (optional) A single value to record as the source of these data points, e.g. "NCBI_browser". By default, the field 'BioProject' will be used to indicate the source for each row in the input file, but if this is missing or you want to override it with a single value for all samples, you may provide a source name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of coalesce if_else mutate relocate rename
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The biological sample identifier (renamed from `BioSample` or specified column).
#' - `spp_pheno`: The species phenotype, formatted using the `as.mo` function.
#' - `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function.
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the `as.mic` function.
#' - `disk`: The disk diffusion measurement (in mm), formatted using the `as.disk` function.
#' - `method`: The AST platform recorded in the input file as the source of the measurement.
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay.
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R), based on the MIC or disk diffusion data.
#' - `pheno_provided`: The original phenotype interpretation provided in the input file.
#' - `source`: The source of each data point (renamed from `BioProject` in the input file, or replaced with a single value passed in as the 'source' parameter).
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
    ast <- ast %>% mutate(method = `Laboratory typing platform`)
  } else {
    cat("Warning: Expected AST platform column 'Laboratory typing platform' not found in input\n")
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
    ast <- ast %>% mutate(pheno_provided = as.sir(`Resistance phenotype`))
  } else {
    cat("Warning: Expected phenotype SIR column 'Resistance phenotype' not found in input\n")
  }

  ast <- interpret_ast(ast, interpret_ecoff = interpret_ecoff, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, species = species, ab = ab)

  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", "pheno_eucast", "pheno_clsi", "ecoff", "guideline", "method", "source", "pheno_provided", "spp_pheno")))

  return(ast)
}


#' Import and Process AST Data from files downloaded from the EBI AMR portal website 
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset that has been downloaded from the EBI AMR portal website (https://www.ebi.ac.uk/amr/data/?view=experiments)
#' Data downloaded from the EBI AMR Portal FTP site (ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/), either directly or via the function `download_ebi()`, is formatted differently and can instead be processed using the `import_ebi_ast_ftp()` function.
#' 
#' This function will process the data, and optionally interpret the results based on MIC or disk diffusion data. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be compressed) and parses relevant columns (antibiotic names, species names, MIC or disk data) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF).
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in EBI antibiogram format. These files can be downloaded from the EBI AMR browser, e.g. https://www.ebi.ac.uk/amr/data/?view=experiments
#' @param sample_col A string indicating the name of the column with sample identifiers. If `NULL`, assume this is 'phenotype-BioSample_ID'.
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the field 'phenotype-organism' will be assumed to specify the species for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the field 'phenotype-antibiotic_name' will be assumed to specify the antibiotic for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @param source (optional) A single value to record as the source of these data points, e.g. "EBI_browser". By default, the field 'phenotype-AMR_associated_publications' will be used to indicate the source for each row in the input file, but if this is missing or you want to override it with a single value for all samples, you may provide a source name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of coalesce if_else mutate relocate rename
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The biological sample identifier (renamed from `phenotype-BioSample_ID` or specified column).
#' - `spp_pheno`: The species phenotype, formatted using the `as.mo` function.
#' - `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function.
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the `as.mic` function.
#' - `disk`: The disk diffusion measurement (in mm), formatted using the `as.disk` function.
#' - `method`: The AST platform recorded in the input file as the source of the measurement.
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay.
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R), based on the MIC or disk diffusion data.
#' - `pheno_provided`: The original phenotype interpretation provided in the input file.
#' - `source`: The source of each data point (renamed from the publications field in the input file, or replaced with a single value passed in as the 'source' parameter).
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
    ast <- ast %>% mutate(method = `phenotype-platform`)
  } else {
    cat("Warning: Expected AST platform column 'phenotype-platform' not found in input\n")
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

  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", "pheno_eucast", "pheno_clsi", "ecoff", "guideline", "method", "source", "pheno_provided", "spp_pheno")))

  return(ast)
}


#' Interpret AST data in a standard format tibble
#'
#' This function applies human EUCAST or CLSI breakpoints, and/or ECOFF, to interpret AST data.
#' @param ast A tibble containing the AST measures in standard AMRgen format, as output by `import_ast`. It must contain assay measurements in columns 'mic' (class mic) and/or 'disk'. Interpretation requires an organism (column 'spp_pheno' of class 'mo', or a single value passed via the 'species' parameter) and an antibiotic (column 'drug_agent' of class 'ab', or a single value passed via the 'ab' parameter).
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the spp_pheno field in the input file will be assumed to specify the species for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the drug_agent field in the input file will be assumed to specify the antibiotic for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr across coalesce mutate
#' @return A copy of the input table, with additional columns:
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R), based on the MIC or disk diffusion data.
#' - `spp_pheno`: The species phenotype, formatted using the `as.mo` function (either taken from the input table, or the single value specified by 'species' parameter).
#' - `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function (either taken from the input table, or the single value specified by 'ab' parameter)..
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
#'   mutate(drug_agent=as.ab(antibiotic)) %>% # convert antibiotic field to 'drug_agent' of class 'ab'
#'   mutate(mic=paste0(sign,MIC)) %>% 
#'   mutate(mic=as.mic(mic)) # create a single 'mic' column of class 'mic'
#'   
#' pheno <- interpret_ast(pheno, species="Escherichia coli")
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
          mutate(across(where(is.mic), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "EUCAST", .names = "pheno_eucast_mic", capped_mic_handling="conservative"))
      }
      if ("disk" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.disk), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "EUCAST", .names = "pheno_eucast_disk"))
      }
      if (("pheno_eucast_mic" %in% colnames(ast)) & ("pheno_eucast_disk" %in% colnames(ast))) {
        ast <- ast %>%
          mutate(pheno_eucast = coalesce(pheno_eucast_mic, pheno_eucast_disk))
      } else if ("pheno_eucast_mic" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_eucast=pheno_eucast_mic)
      } else if ("pheno_eucast_disk" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_eucast=pheno_eucast_disk)
      }
    }
    if (interpret_clsi) {
      if ("mic" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.mic), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "CLSI", .names = "pheno_clsi_mic", capped_mic_handling="conservative"))
      }
      if ("disk" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.disk), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "CLSI", .names = "pheno_clsi_disk"))
      }
      if (("pheno_clsi_mic" %in% colnames(ast)) & ("pheno_clsi_disk" %in% colnames(ast))) {
        ast <- ast %>%
          mutate(pheno_clsi = coalesce(pheno_clsi_mic, pheno_clsi_disk))
      } else if ("pheno_clsi_mic" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_clsi=pheno_clsi_mic)
      } else if ("pheno_eucast_disk" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_clsi_disk=pheno_clsi_disk)
      }
    }
    if (interpret_ecoff) {
      if ("mic" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.mic), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "EUCAST", breakpoint_type = "ECOFF", .names = "ecoff_mic", capped_mic_handling="conservative"))
      }
      if ("disk" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.disk), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "EUCAST", breakpoint_type = "ECOFF", .names = "ecoff_disk"))
      }
      if (("ecoff_mic" %in% colnames(ast)) & ("ecoff_disk" %in% colnames(ast))) {
        ast <- ast %>%
          mutate(ecoff = coalesce(ecoff_mic, ecoff_disk))
      } else if ("ecoff_mic" %in% colnames(ast)) {
        ast <- ast %>% rename(ecoff=ecoff_mic)
      } else if ("ecoff_disk" %in% colnames(ast)) {
        ast <- ast %>% rename(ecoff=ecoff_disk)
      }
    }
  }
  return(ast)
}

#' Import and Process AST Data from an EBI or NCBI antibiogram File
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset in either EBI or NCBI antibiogram format, processes the data, and optionally interprets the results based on MIC or disk diffusion data. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be compressed) and parses relevant columns (antibiotic names, species names, MIC or disk data) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF). If expected columns are not found warnings will be given, and interpretation may not be possible.
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in EBI or NCBI antibiogram format. These files can be downloaded from the EBI AMR web browser (https://www.ebi.ac.uk/amr/data/?view=experiments), EBI FTP site (ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/), or NCBI browser (e.g. https://www.ncbi.nlm.nih.gov/pathogens/ast#Pseudomonas%20aeruginosa), or from EBI using the function `download_ebi()`.
#' @param format A string indicating the format of the data: "ebi" (default), "ebi_web", "ebi_ftp", "ncbi", "vitek", or "whonet". This determines whether the data is passed on to the `import_ebi_ast()` (ebi/ebi_web), `import_ebi_ast_ftp()` (ebi_ftp), `import_ncbi_ast()` (ncbi), `import_vitek_ast()` (vitek), or `import_whonet_ast()` (whonet) function to process.
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the organism field in the input file will be assumed to specify the species for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the antibiotic field in the input file will be assumed to specify the antibiotic for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @param source (optional) A single value to record as the source of these data points, e.g. "EBI_browser". By default, the publications field (for EBI data) or BioProject field (for NCBI data) will be used to indicate the source for each row in the input file, but if this is missing or you want to override it with a single value for all samples, you may provide a source name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of coalesce if_else mutate relocate rename
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The biosample identifier.
#' - `spp_pheno`: The species phenotype, formatted using the `as.mo` function.
#' - `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function.
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the `as.mic` function.
#' - `disk`: The disk diffusion measurement (in mm), formatted using the `as.disk` function.
#' - `method`: The AST platform recorded in the input file as the source of the measurement.
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay.
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R), based on the MIC or disk diffusion data.
#' - `pheno_provided`: The original phenotype interpretation provided in the input file.
#' - `source`: The source of each data point (from the publications or bioproject field in the input file, or replaced with a single value passed in as the 'source' parameter).
#' @export
#' @examples
#' # small example E. coli AST data from NCBI
#' ecoli_ast_raw
#'
#' # import without re-interpreting resistance
#' pheno <- import_ast(ecoli_ast_raw, format="ncbi")
#' head(pheno)
#'
#' # import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
#' pheno <- import_ast(ecoli_ast_raw, format="ncbi", interpret_eucast = TRUE, interpret_ecoff = TRUE)
#' head(pheno)
import_ast <- function(input, format = "ebi", interpret_eucast = FALSE,
                       interpret_clsi = FALSE, interpret_ecoff = FALSE,
                       species = NULL, ab = NULL, source=NULL) {
  if (format %in% c("ebi_web", "ebi")) {
    cat("Reading in as EBI AST format downloaded from the web portal\n")
    ast <- import_ebi_ast(input, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, interpret_ecoff = interpret_ecoff, species = species, ab = ab)
  }
  
  if (format == "ebi_ftp") {
    cat("Reading in as EBI AST format downloaded from the FTP portal\n")
    ast <- import_ebi_ast_ftp(input, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, interpret_ecoff = interpret_ecoff)
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
                            species = species, ab = ab, source = source)
  }

  if (format == "whonet") {
    cat("Reading in as WHONET AST format\n")
    ast <- import_whonet_ast(input,
                             interpret_eucast = interpret_eucast,
                             interpret_clsi = interpret_clsi,
                             interpret_ecoff = interpret_ecoff,
                             species = species, ab = ab, source = source)
  }

  if (!is.null(source)) { ast <- ast %>% mutate(source=source)}
  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", "pheno_eucast", "pheno_clsi", "ecoff", "guideline", "method", "source", "pheno_provided", "spp_pheno")))

  return(ast)
}


#' Import and Process AST Data from a generic format
#'
#' This function attempts to import antibiotic susceptibility testing (AST) data in long-form antibiogram format (one row per sample and test), suitable for downstream use with AMRgen analysis functions. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be compressed) and parses relevant columns (antibiotic names, species names, MIC or disk data, S/I/R calls) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF). If expected columns are not found warnings will be given, and interpretation may not be possible.
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in long-form antibiogram format (one row per sample and test). This might be a file containing the content of an EBI/NCBI AST dataset previously processed using `import_ebi_ast()`, `import_ncbi_ast()`, or `import_ast()` functions, or files with a similar format/structure but with different column names.
#' @param sample_col (optional, default 'id') Name of the input data column that provides the sample name. If the 'rename' parameter is set to TRUE, this column will be renamed as 'id'.
#' @param species_col (optional, default 'species') Name of the input data column that provides a species name. If provided, this column will be converted to micro-organism class 'mo' via `as.mo()`. If the 'rename' parameter is set to TRUE, this column will also be renamed as 'spp_pheno'. If interpretation is switched on, this column will be used to identify the appropriate  breakpoints for interpretation of each row in the data table.
#' @param species (optional) Name of the single species to which all samples belong. Use this if you want to interpret assay measurements but the input file does not contain a column indicating the species for each sample (called 'species' or a name specified by the `species_col` parameter).
#' @param ab_col (optional, default 'drug_agent') Name of the input data column that provides a drug name. If provided, this column will be converted to antibiotic class 'ab' via `as.ab()`. If the 'rename' parameter is set to TRUE, this column will also be renamed as 'drug_agent'. If interpretation is switched on, this column will be used to identify the appropriate breakpoints for interpretation of each row in the data table.
#' @param ab (optional) Name of a single antibiotic to use for phenotype interpretation. Use this if you want to interpret assay measurements but the input file does not contain a column indicating the drug for each sample (called 'drug_agent' or a name specified by the `ab_col` parameter).
#' @param mic_col (optional, default 'mic') Name of the input data column that provides MIC measurements. If provided, this column will be converted to MIC class 'mic' via `as.mic()`. If the 'rename' parameter is set to TRUE, this column will also be renamed as 'mic'. If interpretation is switched on, the MIC values will be interpreted against clinical breakpoints.
#' @param disk_col (optional, default 'disk') Name of the input data column that provides disk diffusion zone measurements. If provided, this column will be converted to disk diffusion class 'disk' via `as.disk()`. If the 'rename' parameter is set to TRUE, this column will also be renamed as 'disk'. If interpretation is switched on, the zone values will be interpreted against clinical breakpoints.
#' @param pheno_cols (optional, default `c("ecoff", "pheno_eucast", "pheno_clsi", "pheno_provided")`) Name of the input data column/s that provides disk diffusion zone measurements (as a character vector, or single string for a single column). If provided, these columns will be converted to SIR class 'sir' via `as.sir()`.
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @param rename_cols A logical value (default is TRUE). If `TRUE`, the function will rename the provided columns (specified by `ab_col`, `mic_col`, `disk_col`, `species_col`, `id_col`) to the default names expected by AMRgen functions ('drug_agent', 'mic', 'disk', 'spp_pheno', 'id'), to match those output by the other `import_ast()` functions.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of mutate relocate
#' @importFrom rlang is_string :=
#' @return A data frame with the processed AST data, including additional columns:
#' @export
#' @examples
#' \dontrun{
#' # import and process AST data from EBI, write formatted data to file for later use
#' pheno <- import_ebi_ast("EBI_AMR_data.csv.gz")
#' write_tsv(pheno, file="EBI_AMR_data_processed.tsv.gz", 
#'             interpret_eucast = TRUE, interpret_ecoff = TRUE)
#'
#' # read stored data and format the columns to the correct classes
#' pheno <- format_ast("EBI_AMR_data_processed.tsv.gz")
#'
#' # read in unprocessed E. coli AST data from non-standard format and interpret
#' pheno <- format_ast("AMR_data.tsv", sample_col="STRAIN", species="E. coli",
#' 						ab_col="Antibiotic", mic_col="MIC (mg/L)", 
#' 						interpret_eucast=TRUE, interpret_ecoff=TRUE)
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
                       interpret_eucast = FALSE, 
                       interpret_clsi = FALSE, 
                       interpret_ecoff = FALSE,
                       rename_cols = TRUE) {
  
  ast <- process_input(input)
  
  if (!is.null(species)) {
    cat(paste("Adding new micro-organism column 'spp_pheno' (class 'mo') with constant value",species,"\n"))
    ast$spp_pheno <- as.mo(species)
  }
  
  if (!is.null(species_col)) {
    if (species_col %in% colnames(ast)) {
      ast <- ast %>% mutate(!!sym(species_col):=as.mo(!!sym(species_col)))
      cat(paste("Parsing column", species_col,"as micro-organism (class 'mo')\n"))
      if (rename_cols & species_col!="spp_pheno") {
        ast <- ast %>% rename(spp_pheno=!!sym(species_col))
        cat(paste("Renaming column", species_col,"to standard name 'spp_pheno'\n"))
      } 
      if ((interpret_ecoff | interpret_eucast | interpret_clsi) & !("spp_pheno" %in% colnames(ast))) {
        ast <- ast %>% mutate(spp_pheno=!!sym(species_col)) # we need a column named spp_pheno for interpretation
      }
    }
    else { cat(paste("Could not find species_col", species_col, "in input table")) }
  }
  
  if (!is.null(ab)) {
    cat(paste("Adding new antibiotic column 'drug_agent' (class 'ab') with constant value",ab,"\n"))
    ast$drug_agent <- as.ab(ab)
  }
  
  if (!is.null(ab_col)) {
    if (ab_col %in% colnames(ast)) {
      ast <- ast %>% mutate(!!sym(ab_col):=as.ab(!!sym(ab_col)))
      cat(paste("Parsing column", ab_col,"as antibiotic (class 'ab')\n"))
      if (rename_cols & ab_col!="drug_agent") {
        ast <- ast %>% rename(drug_agent=!!sym(ab_col))
        cat(paste("Renaming column", ab_col,"to standard name 'drug_agent'\n"))
      }
    }
    else { cat(paste("Could not find ab_col", ab_col, "in input table")) }
  }
  
  if (!is.null(mic_col)) {
    if (mic_col %in% colnames(ast)) {
      ast <- ast %>% mutate(!!sym(mic_col):=as.mic(!!sym(mic_col)))
      cat(paste("Parsing column", mic_col,"as class 'mic'\n"))
      if (rename_cols & mic_col!="mic") {
        ast <- ast %>% rename(mic=!!sym(mic_col))
        cat(paste("Renaming column", mic_col,"to standard name 'mic'\n"))
      }
    }
    else { cat(paste("Could not find mic_col", mic_col, "in input table")) }
  }
  
  if (!is.null(disk_col)) {
    if (disk_col %in% colnames(ast)) {
      ast <- ast %>% mutate(!!sym(disk_col):=as.disk(!!sym(disk_col)))
      cat(paste("Parsing column", disk_col,"as class 'disk'\n"))
      if (rename_cols & disk_col!="disk") {
        ast <- ast %>% rename(disk=!!sym(disk_col))
        cat(paste("Renaming column", disk_col,"to standard name 'disk'\n"))
      }
    }
    else { cat(paste("Could not find disk_col", disk_col, "in input table")) }
  }
  
  if (!is.null(pheno_cols)) {
    if (is_string(pheno_cols)) { # single column
      if (pheno_cols %in% colnames(ast)) {
        ast <- ast %>% mutate(!!sym(pheno_cols):=as.sir(!!sym(pheno_cols)))
        cat(paste("Parsing column", pheno_cols,"as class 'sir'\n"))
      }
      else { cat(paste("Could not find pheno_cols", pheno_cols, "in input table")) }
    }
    else {
      for (pheno_col in pheno_cols) {
        if (!is.null(pheno_col)) {
          if (pheno_col %in% colnames(ast)) {
            ast <- ast %>% mutate(!!sym(pheno_col):=as.sir(!!sym(pheno_col)))
            cat(paste("Parsing column", pheno_col,"as class 'sir'\n"))
          }
          else { cat(paste("Could not find pheno_col", pheno_col, "in input table")) }
        }
      }
    }
  }
  
  ast <- interpret_ast(ast, interpret_ecoff = interpret_ecoff, 
                       interpret_eucast = interpret_eucast, 
                       interpret_clsi = interpret_clsi, 
                       species = species, # note column spp_pheno will be used otherwise
                       ab = ab)
  
  if (rename_cols) {ast <- ast %>% rename(id=!!sym(sample_col))}
  
  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", ab_col, mic_col, disk_col, pheno_cols, 
                                   "pheno_eucast", "pheno_clsi", "ecoff", "guideline", 
                                   "method", "source", "pheno_provided", "spp_pheno", species_col
  )))
  
  return(ast)
}


#' Import and Process AST Data files retrieved from the EBI AMR portal FTP site
#'
#' This function will import antibiotic susceptibility testing (AST) data suitable for downstream use with AMRgen analysis functions. The expected input is phenotype data retrieved from the EBI AMR Portal FTP site (ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) either directly or via the function `download_ebi()`.
#' Note that files downloaded from the EBI AMR Portal web browser (https://www.ebi.ac.uk/amr/data/?view=experiments) are formatted differently and can be imported with the function `import_ebi_ast()`
#' @param input A string representing the input dataframe, or a path to an input file, to be processed.
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of mutate relocate
#' @importFrom rlang is_string :=
#' @return A data frame with the processed AST data, including additional columns:
#' @export
#' @examples
#' \dontrun{
#' # download Salmonella phenotype data from EBI
#' pheno_salmonella <- download_ebi(genus="Salmonella")
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
    mutate(mic=if_else(measurement_units=="mg/L", 
                       paste0(measurement_sign, measurement),
                       NA)) %>%
    mutate(disk=if_else(laboratory_typing_method=="disk diffusion", 
                        paste0(measurement_sign, measurement),
                        NA)) %>%
    mutate(pheno_provided=if_else(resistance_phenotype=="intermediate", 
                                  "I", 
                                  resistance_phenotype)) %>%
    rename(method=platform, source=AMR_associated_publications, guideline=ast_standard) %>%
    mutate(source=if_else(is.na(source), "EBI", paste0("EBI:",source))) %>%
    format_ast(sample_col="BioSample_ID", 
               species_col="species", 
               ab_col="antibiotic_name", 
               pheno_cols="pheno_provided",
               interpret_eucast = interpret_eucast,
               interpret_clsi = interpret_clsi,
               interpret_ecoff = interpret_ecoff)

  return(ast)
}


#' Import and Process AST Data from Vitek Output Files
#'
#' This function imports AST data from Vitek instrument output files (wide CSV format)
#' and converts it to the standardized long-format used by AMRgen.
#'
#' @param input A dataframe or path to a CSV file containing Vitek AST output data
#' @param sample_col Column name for sample identifiers. Default: "Lab ID"
#' @param source Optional source value to record for all data points
#' @param species Optional species override for phenotype interpretation
#' @param ab Optional antibiotic override for phenotype interpretation
#' @param use_expertized Use expertized SIR (TRUE, default) or instrument SIR (FALSE)
#' @param interpret_eucast Interpret against EUCAST breakpoints
#' @param interpret_clsi Interpret against CLSI breakpoints
#' @param interpret_ecoff Interpret against ECOFF values
#' @param include_dates Include collection_date and testing_date in output
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr across any_of case_when coalesce mutate relocate rename
#' @importFrom tidyr pivot_longer matches
#' @importFrom stringr str_match
#' @importFrom rlang sym
#' @return Standardized AST data frame
#' @export
import_vitek_ast <- function(input,
                             sample_col = "Lab ID",
                             source = NULL,
                             species = NULL,
                             ab = NULL,
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

  # Identify antibiotic MIC columns (CODE-Name pattern, excluding -Other-*)
  all_cols <- colnames(ast)
  mic_cols <- all_cols[grepl("^[A-Z0-9]+-", all_cols) &
                        !grepl("-Other-Instrument$|-Other-Expertized$", all_cols)]

  if (length(mic_cols) == 0) {
    stop("No antibiotic columns found in Vitek format")
  }

  # Extract antibiotic codes and names from column headers
  ab_codes <- stringr::str_match(mic_cols, "^([A-Z0-9]+)-")[,2]
  ab_names <- stringr::str_match(mic_cols, "^[A-Z0-9]+-(.*)$")[,2]

  # Metadata columns to preserve
  metadata_cols <- c("Lab ID", "Isolate Number", "Patient Name", "Patient ID",
                     "Patient Location", "Specimen Type", "Specimen Source",
                     "Collection Date", "Testing Date", "Organism Name",
                     "Organism Code", "Bio Number", "Percent Probability",
                     "ID Confidence", "Selected BP Infection Site")
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
    mutate(across(any_of(all_ab_cols), as.character))

  # Pivot to long format
  ast_long <- ast %>%
    tidyr::pivot_longer(
      cols = tidyr::matches("^[A-Z0-9]+_(mic|sir_inst|sir_exp)$"),
      names_to = c("ab_code", ".value"),
      names_pattern = "^([A-Z0-9]+)_(mic|sir_inst|sir_exp)$"
    )

  # Add antibiotic full name and parse values
  ast_long <- ast_long %>%
    mutate(ab_name = ab_lookup[ab_code]) %>%
    mutate(mic = case_when(
      is.na(mic) | mic == "" | mic == "-" | mic == "NEG" ~ NA_character_,
      TRUE ~ mic
    )) %>%
    mutate(mic = as.mic(mic)) %>%
    mutate(disk = as.disk(NA))  # Vitek is MIC-based only

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
  ast_long <- ast_long %>% mutate(pheno_provided = as.sir(pheno_provided))

  # Parse organism and antibiotic
  if ("Organism Name" %in% colnames(ast_long)) {
    ast_long <- ast_long %>% mutate(spp_pheno = as.mo(`Organism Name`))
  }

  # Try code first, then name for drug_agent
  ast_long <- ast_long %>%
    mutate(drug_agent_code = as.ab(ab_code)) %>%
    mutate(drug_agent_name = as.ab(ab_name)) %>%
    mutate(drug_agent = coalesce(drug_agent_code, drug_agent_name))

  # Rename sample column and add standard columns
  ast_long <- ast_long %>%
    rename(id = !!sym(sample_col)) %>%
    mutate(method = "VITEK", guideline = NA_character_)

  # Add source
  if (!is.null(source)) {
    ast_long <- ast_long %>% mutate(source = source)
  } else {
    ast_long <- ast_long %>% mutate(source = "VITEK")
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
                            species = species, ab = ab)

  # Reorder columns
  ast_long <- ast_long %>%
    relocate(any_of(c("id", "drug_agent", "mic", "disk",
                      "pheno_eucast", "pheno_clsi", "ecoff",
                      "guideline", "method", "source",
                      "pheno_provided", "spp_pheno")))

  return(ast_long)
}


#' Import and Process AST Data from WHONET Output Files
#'
#' This function imports AST data from WHONET software output files (wide CSV format)
#' and converts it to the standardised long-format used by AMRgen.
#'
#' @param input A dataframe or path to a CSV file containing WHONET AST output data
#' @param sample_col Column name for sample identifiers. Default: "Identification number"
#' @param source Optional source value to record for all data points
#' @param species Optional species override for phenotype interpretation
#' @param ab Optional antibiotic override for phenotype interpretation
#' @param interpret_eucast Interpret against EUCAST breakpoints
#' @param interpret_clsi Interpret against CLSI breakpoints
#' @param interpret_ecoff Interpret against ECOFF values
#' @param include_patient_info Include patient demographic columns in output
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr across any_of case_when coalesce mutate relocate rename
#' @importFrom tidyr pivot_longer matches
#' @importFrom stringr str_match
#' @importFrom rlang sym
#' @return Standardised AST data frame
#' @export
import_whonet_ast <- function(input,
                              sample_col = "Identification number",
                              source = NULL,
                              species = NULL,
                              ab = NULL,
                              interpret_eucast = FALSE,
                              interpret_clsi = FALSE,
                              interpret_ecoff = FALSE,
                              include_patient_info = FALSE) {

  ast <- process_input(input)

  # Validate sample column exists
  if (!(sample_col %in% colnames(ast))) {
    stop(paste("Invalid column name:", sample_col))
  }

  # Identify antibiotic columns (CODE_METHOD pattern, e.g., AMP_ND10, CIP_ED5, AMP_EE)
  all_cols <- colnames(ast)
  ab_cols <- all_cols[grepl("^[A-Z]{2,4}_[A-Z]{2}", all_cols)]

  if (length(ab_cols) == 0) {
    stop("No antibiotic columns found in WHONET format")
  }

  # Extract antibiotic codes from column headers
  ab_codes <- stringr::str_match(ab_cols, "^([A-Z]{2,4})_")[,2]

  # Metadata columns to preserve
  metadata_cols <- c("Identification number", "Specimen number", "Organism",
                     "Country", "Laboratory", "Specimen date", "Specimen type",
                     "Isolate number", "Organism type", "Date of data entry")
  if (include_patient_info) {
    metadata_cols <- c(metadata_cols, "Last name", "First name", "Sex",
                       "Age", "Age category", "Date of admission")
  }
  metadata_cols <- metadata_cols[metadata_cols %in% all_cols]

  # Convert antibiotic columns to character to avoid type conflicts
  ast <- ast %>%
    mutate(across(any_of(ab_cols), as.character))

  # Pivot to long format
  ast_long <- ast %>%
    tidyr::pivot_longer(
      cols = any_of(ab_cols),
      names_to = "ab_col",
      values_to = "sir_value"
    )

  # Parse antibiotic code and method from column name
  ast_long <- ast_long %>%
    mutate(ab_code = stringr::str_match(ab_col, "^([A-Z]{2,4})_")[,2]) %>%
    mutate(test_method = stringr::str_match(ab_col, "^[A-Z]{2,4}_(.*)$")[,2])

  # Parse SIR values
  ast_long <- ast_long %>%
    mutate(pheno_provided = case_when(
      sir_value %in% c("S", "I", "R") ~ sir_value,
      TRUE ~ NA_character_
    )) %>%
    mutate(pheno_provided = as.sir(pheno_provided))

  # Set MIC and disk to NA (WHONET provides interpretations, not raw values)
  ast_long <- ast_long %>%
    mutate(mic = as.mic(NA)) %>%
    mutate(disk = as.disk(NA))

  # Parse organism
  if ("Organism" %in% colnames(ast_long)) {
    ast_long <- ast_long %>% mutate(spp_pheno = as.mo(Organism))
  }

  # Parse antibiotic code to drug_agent
  ast_long <- ast_long %>%
    mutate(drug_agent = as.ab(ab_code))

  # Rename sample column and add standard columns
  ast_long <- ast_long %>%
    rename(id = !!sym(sample_col)) %>%
    mutate(method = "WHONET", guideline = NA_character_)

  # Add source
  if (!is.null(source)) {
    ast_long <- ast_long %>% mutate(source = source)
  } else {
    ast_long <- ast_long %>% mutate(source = "WHONET")
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
                            species = species, ab = ab)

  # Reorder columns
  ast_long <- ast_long %>%
    relocate(any_of(c("id", "drug_agent", "mic", "disk",
                      "pheno_eucast", "pheno_clsi", "ecoff",
                      "guideline", "method", "source",
                      "pheno_provided", "spp_pheno")))

  return(ast_long)
}
