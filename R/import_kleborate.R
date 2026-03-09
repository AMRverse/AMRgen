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

#' Import and Process Kleborate Results
#'
#' This function imports and processes Kleborate (https://github.com/klebgenomics/Kleborate) results, extracting antimicrobial resistance (AMR) elements and mapping them to standardised drug classes. 
#' @param input_table A character string specifying a dataframe or path to the Kleborate results table (TSV format).
#' @param sample_col A character string specifying the column that identifies samples in the dataset (default `strain`).
#' @importFrom dplyr filter left_join mutate select bind_rows rename_with
#' @importFrom tidyr separate_longer_delim separate
#' @importFrom rlang sym
#' @return A tibble containing the processed AMR determinants and drug classes that is AMRgen compatible. 
#' @details
#' The function performs the following steps:
#' - Reads the Kleborate output table.
#' - Transforms Kleborate output into long form (i.e., one AMR determinant per row).
#' - Maps Kleborate drug classes to standardised drug class.
#' This processing ensures compatibility with downstream AMRgen analysis workflows.
#' @export
#' @examples
#' \dontrun{
#' # small example Kleborate data
#' data(kleborate_raw)
#' kleborate_raw
#'
#' # import first few rows of this data frame and parse it as AMRfp data
#' geno <- import_kleborate(kleborate_raw %>% head(n = 10), "strain")
#' geno
#' }


import_kleborate <- function(input_table,
                             sample_col = "strain") {
  in_table <- process_input(input_table)
  
  non_mutation_AMR_table <- make_non_mutation_AMR_table(in_table, c("AGly_acquired",
                                                                    "Col_acquired",
                                                                    "Fcyn_acquired",
                                                                    "Flq_acquired",
                                                                    "Gly_acquired",
                                                                    "MLS_acquired",
                                                                    "Phe_acquired",
                                                                    "Rif_acquired",
                                                                    "Sul_acquired",
                                                                    "Tet_acquired",
                                                                    "Tgc_acquired",
                                                                    "Tmt_acquired",
                                                                    "Bla_acquired",
                                                                    "Bla_inhR_acquired",
                                                                    "Bla_ESBL_acquired",
                                                                    "Bla_ESBL_inhR_acquired",
                                                                    "Bla_Carb_acquired",
                                                                    "Bla_chr"))
  
  mutation_AMR_table <- make_mutation_AMR_table(in_table, c("Omp_mutations", 
                                                            "Col_mutations", 
                                                            "Flq_mutations"))
  
  AMR_table <- dplyr::bind_rows(non_mutation_AMR_table, mutation_AMR_table)
  AMR_table$drug_agent <- NA
  
  if (!"species" %in% colnames(in_table)) {
    stop("Input file lacks the expected column: species")}
  else {
    species_table <- in_table %>% 
      select(strain, species)
    
    species_table$genus <- strsplit(species_table$species, " ")[[1]][1]
    species_table$organism <- species_table$species
  }
  
  kleborate_AMR_table <- dplyr::left_join(AMR_table, species_table, by = "strain")
  return(kleborate_AMR_table)
}


# Function to transform AMR gene (non-mutation) columns into long form
make_non_mutation_AMR_table <- function(input_table, column_names) {
  
  in_table <- process_input(input_table)
  
  results <- list()
  
  for (column_name in column_names) {
    
    if (!column_name %in% colnames(in_table)) {
      stop(paste("Input file lacks the expected column:", column_name))
    }
    
    if (!any(in_table[[column_name]] != "-" & !is.na(in_table[[column_name]]))) {
      message(paste("No", column_name, "in any genomes"))
      next
    }
    
    df <- in_table %>% 
      dplyr::select(strain, !!rlang::sym(column_name)) %>% 
      dplyr::filter(.data[[column_name]] != "-") %>%
      tidyr::separate_longer_delim(
        cols = !!rlang::sym(column_name),
        delim = ";"
      ) %>%
      dplyr::rename(gene = !!rlang::sym(column_name)) %>%
      # remove all '^' from gene which indicates exact protein sequence match, inexact nucleotide match 
      dplyr::mutate(gene = str_remove_all(gene, "\\^"))
    
    df$marker <- df$gene
    df$marker.label <- df$gene
    df$mutation <- "-"
    df$element_type <- "AMR"
    df$element_subtype <- "AMR"
    
    # Define drug class classifications
    if (grepl("^AGly", column_name)) {
      df$drug_class <- "Aminoglycosides"
      df$class <- "AMINOGLYCOSIDE"
      df$subclass <- NA
      
    } else if (grepl("^Col_acquired", column_name)) {
      df$drug_class <- "Polymyxins"
      df$class <- "COLISTIN"
      df$subclass <- "COLISTIN"
      
    } else if (grepl("^Fcyn_acquired", column_name)) {
      df$drug_class <- "Phosphonics"
      df$class <- "FOSFOMYCIN"
      df$subclass <- "FOSFOMYCIN"
      
    } else if (grepl("^Flq_acquired", column_name)) {
      df$drug_class <- "Quinolones"
      df$class <- "QUINOLONE"
      df$subclass <- "QUINOLONE"
      
    } else if (grepl("^Gly_acquired", column_name)) {
      df$drug_class <- "Glycopeptides"
      df$class <- NA
      df$subclass <- NA
      
    } else if (grepl("^MLS_acquired", column_name)) {
      df$drug_class <- "Macrolides/lincosamides"
      df$class <- "MACROLIDE"
      df$subclass <- NA
      
    } else if (grepl("^Phe_acquired", column_name)) {
      df$drug_class <- "Phenicols"
      df$class <- "PHENICOL"
      df$subclass <- NA
      
    } else if (grepl("^Rif_acquired", column_name)) {
      df$drug_class <- "Antimycobacterials"
      df$class <- "RIFAMYCIN"
      df$subclass <- "RIFAMYCIN"
      
    } else if (grepl("^Sul_acquired", column_name)) {
      df$drug_class <- "Other antibacterials"
      df$class <- "SULFONAMIDE"
      df$subclass <- "SULFONAMIDE"
      
    } else if (grepl("^Tet_acquired", column_name)) {
      df$drug_class <- "Tetracyclines"
      df$class <- "TETRACYCLINE"
      df$subclass <- "TETRACYCLINE"
      
    } else if (grepl("^Tgc_acquired", column_name)) {
      df$drug_class <- "Tetracyclines"
      df$class <- "TETRACYCLINE"
      df$subclass <- "TIGECYCLINE"
      
    } else if (grepl("^Tmt_acquired", column_name)) {
      df$drug_class <- "Trimethoprims"
      df$class <- "TRIMETHOPRIM"
      df$subclass <- "TRIMETHOPRIM"
      
    } else if (grepl("^Bla_acquired", column_name)) {
      df$drug_class <- "Beta-lactams"
      df$class <- "BETA-LACTAM"
      df$subclass <- "BETA-LACTAM"
      
    } else if (grepl("^Bla_inhR_acquired", column_name)) {
      df$drug_class <- "Beta-lactams/penicillins"
      df$class <- "BETA-LACTAM"
      df$subclass <- "BETA-LACTAM"
      
    } else if (grepl("^Bla_ESBL_acquired", column_name)) {
      df$drug_class <- "Cephalosporins (3rd gen.)"
      df$class <- "BETA-LACTAM"
      df$subclass <- "CEPHALOSPORIN"
      
    } else if (grepl("^Bla_ESBL_inhR_acquired", column_name)) {
      df$drug_class <- "Cephalosporins (3rd gen.)"
      df$class <- "BETA-LACTAM"
      df$subclass <- "CEPHALOSPORIN"
      
    } else if (grepl("^Bla_Carb_acquired", column_name)) {
      df$drug_class <- "Carbapenems"
      df$class <- "BETA-LACTAM"
      df$subclass <- "CARBAPENEM"
      
    } else if (grepl("^Bla_chr", column_name)) {
      df$drug_class <- "Beta-lactams/penicillins"
      df$class <- "BETA-LACTAM"
      df$subclass <- "BETA-LACTAM"
    }
    
    results[[column_name]] <- df
  }
  
  if (length(results) == 0) {
    return(NULL)
  }
  
  dplyr::bind_rows(results)
}

# Function to transform AMR mutation columns into long form
make_mutation_AMR_table <- function(input_table, column_names) {
  
  in_table <- process_input(input_table)
  
  results <- list()
  
  for (column_name in column_names) {
    
    if (!column_name %in% colnames(in_table)) {
      stop(paste("Input file lacks the expected column:", column_name))
    }
    
    if (!any(in_table[[column_name]] != "-" & !is.na(in_table[[column_name]]))) {
      message(paste("No", column_name, "in any genomes"))
      next
    }
    
    df <- in_table %>%
      dplyr::select(strain, !!rlang::sym(column_name)) %>%
      dplyr::filter(.data[[column_name]] != "-") %>%
      tidyr::separate_longer_delim(
        cols = !!rlang::sym(column_name),
        delim = ";"
      ) %>%
      dplyr::rename(marker.label = !!rlang::sym(column_name)) %>%
      # To split protein/nucleotide mutations separated by ':' and trunctions separated by '_' 
      tidyr::separate(
        marker.label,
        into = c("gene", "mutation"),
        sep = "[:_]",
        remove = FALSE,
        extra = "merge",
        fill = "right"
      ) %>%
      # For pmrB$ or mgrB$ which indicate a mutation in the start codon that may disrupt translation 
      dplyr::mutate(
        gene = ifelse(stringr::str_ends(marker.label, "\\$"),
                      stringr::str_remove(marker.label, "\\$"),
                      gene),
        mutation = ifelse(stringr::str_ends(marker.label, "\\$"),
                          "$",
                          mutation)
      )
    
    df$marker <- paste(df$gene, df$mutation, sep = "_")
    df$element_type <- "AMR"
    df$element_subtype <- "POINT"
    
    # Define drug class classifications
    if (grepl("^Omp_mutations", column_name)) {
      df$drug_class <- "Carbapenems"
      df$class <- "BETA-LACTAM"
      df$subclass <- "CARBAPENEM"
      
    } else if (grepl("^Col_mutations", column_name)) {
      df$drug_class <- "Polymyxins"
      df$class <- "COLISTIN"
      df$subclass <- "COLISTIN"
      
    } else if (grepl("^Flq_mutations", column_name)) {
      df$drug_class <- "Quinolones"
      df$class <- "QUINOLONE"
      df$subclass <- "QUINOLONE"
    }
    
    results[[column_name]] <- df
  }
  
  if (length(results) == 0) {
    return(NULL)
  }
  
  dplyr::bind_rows(results)
}

