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

#' Import and Process Resistance Gene Identifier (RGI) Results
#'
#' This function imports and processes genotyping results from the Resistance Gene Identifier (RGI, https://github.com/arpcard/rgi), extracting antimicrobial resistance determinants and mapping them to standardised drug classes/antibiotics.
#' @param input_table A character string specifying a dataframe or path to the RGI results table (TSV format).
#' @param orf_id_col A character string specifying the column that identifies open reading frame ID (ORF_ID) in the dataset (default `ORF_ID`). This column includes the sample ID in addition to contig / genomic location information about where the resistance determinant was identified and is a default output of RGI.  
#' @param sample_id_sep A character string specifying the separater by which the sample ID is separated from open reading frame ID (ORF_ID) details. For example: in the `ORF_ID` column, "SAMEA3498968.fasta.txt:1_96 # 109511 # 110635....", the sample ID separater is `.fasta.txt:` where SAMEA3498968 is the Sample ID and the remaining text `1_96 # 109511 # 110635...` is the open reading frame information. Default is `.fasta.txt:`.  
#' @param model_col A character string specifying the column that identifies model type identified by RGI (default `Model_type`).
#' @param antibiotic_col Character string specifying the antibiotic column (default `Antibiotic`).
#' @param class_col Character string specifying the drug class column (default `Drug Class`).
#' @param exclude_loose Logical indicating whether to exclude Loose hits (AMR markers that fall below a curated bitscore cutoff as defined by CARD/RGI). Default `TRUE`, which excludes Loose hits. 
#' @param rgi_short_name A tibble containing a reference table mapping model IDs (from CARD/RGI) to shortened model names as provided by CARD (https://card.mcmaster.ca/download in aro_index.tsv). Defaults to `rgi_short_name_table`, which is provided internally.
#' @param rgi_drugs A tibble containing a reference table mapping CARD drug class / drug agents to standardised drug classes/names. Defaults to `rgi_drugs_table`, which is provided internally.
#' @importFrom AMR as.ab ab_group
#' @importFrom dplyr filter left_join mutate select bind_rows
#' @importFrom tidyr separate_longer_delim
#' @importFrom stringr str_remove_all
#' @return A tibble containing the processed AMR determinants and drug classes that is AMRgen compatible. The output retains the original columns from the RGI output along with the newly mapped variables.
#' @details
#' The function performs the following steps:
#' - Reads the RGI output table.
#' - Transforms RGI output into long form (i.e., one AMR determinant AND drug class / antibiotic per row).
#' - Maps CARD drug classes and antibiotics to standardised names.
#' This processing ensures compatibility with downstream AMRgen analysis workflows.
#' @export
#' @examples
#' # example RGI data (including Perfect, Strict, and Loose hits)
#' rgi_raw
#'
#' # import and parse it as AMRfp data using sample_id_sep=`_genomic.fna.txt:` and include Loose hits
#' import_rgi(rgi_raw, sample_id_sep="_genomic.fna.txt:", exclude_loose = FALSE)
#' 
#' # example RGI data from EuSCAPE project (including only Perfect and Strict hits)
#' rgi_EuSCAPE_raw
#' 
#' # import and parse it as AMRfp data using defaults sample_id_sep=`.fasta.txt:` and exclude Loose hits (e.g., exclude_loose = `TRUE`)
#' import_rgi(rgi_EuSCAPE_raw)


import_rgi <- function(input_table,
                         orf_id_col = "ORF_ID",
                         sample_id_sep = ".fasta.txt:",
                         model_col = "Model_type",
                         antibiotic_col = "Antibiotic",
                         class_col = "Drug Class",
                         exclude_loose = TRUE,
                         rgi_short_name = rgi_short_name_table,
                         rgi_drugs = rgi_drugs_table) {
  in_table <- process_input(input_table)

  in_table <- left_join(in_table, rgi_short_name, by=c("Model_ID"="Model ID")) %>%
    mutate(id = sub(paste0(sample_id_sep, ".*"), "", .data[[orf_id_col]]))
  
  in_table[(in_table == "n/a")|(in_table == "")] <- NA
  
  if(exclude_loose) { # exclude Loose hits
    geno_table <- in_table %>% filter(`Cut_Off`!="Loose")
  }
  else{ # include all hits
    geno_table <- in_table 
  }
  
  # AMR marker name assignment
  if ("CARD.Short.Name" %in% colnames(geno_table)) { # Use CARD Short Name for marker label
    geno_table <- geno_table %>% mutate(marker = `CARD Short Name`)
  } else if ("Best_Hit_ARO" %in% colnames(geno_table)) {# Use Best_Hit_ARO for marker label (noting that the names could be very long)
    geno_table <- geno_table %>% mutate(marker = `Best_Hit_ARO`)
  } else {
    stop(paste("Input file lacks the expected column: `CARD Short Name` OR Best_Hit_ARO \n"))
  }
  
  # Assign variation type based on RGI 'Model_type' column
  if (model_col %in% colnames(geno_table)) {
    geno_table <- geno_table %>%
      mutate(`variation type` = case_when(
        !!sym(model_col) == "protein homolog model" ~ "Gene presence detected",
        !!sym(model_col) %in% c("protein variant model", "protein overexpression model") ~ "Protein variant detected",
        !!sym(model_col) == "rRNA gene variant model" ~ "Nucleotide variant detected",
        TRUE ~ NA
      ))
    } 
  else {
    cat("Need method column: Model_type", "\n")
  }
  
  # Check mutation column exists and reformat table to long form - one mutation per row
  if ("SNPs_in_Best_Hit_ARO" %in% colnames(geno_table)) { 
    geno_table <- geno_table %>% 
      rename(mutation = SNPs_in_Best_Hit_ARO)  %>% 
      separate_longer_delim(mutation, delim = ", ")

  } else {
    stop(paste("Input file lacks the expected column: SNPs_in_Best_Hit_ARO \n"))
  }
  
  # create AMRrules style label with gene:mutation
    variant_models <- c(
      "protein variant model",
      "protein overexpression model",
      "rRNA gene variant model")
    
    geno_table_label <- geno_table %>%
      mutate(
        marker.label = case_when(
          !!sym(model_col) %in% variant_models ~ if_else(
            !is.na(mutation),
            paste(`CARD Short Name`, mutation, sep = "_"),
            paste0(`CARD Short Name`, ":-")
          ),
          !!sym(model_col) == "protein homolog model" ~ if_else(
            Cut_Off == "Perfect",
            `CARD Short Name`,
            paste0(`CARD Short Name`, ":-") # Strict or Loose hit           
            ),
          TRUE ~ `CARD Short Name`
        )
      )
    
    # Identify any drug classes / antibiotics we _know_ aren't in the AMR package, using the internal data
    # Join introduces these as new drug_internal or drug_class_internal column
    # Standardizing antibiotic names & drug class first, if no antibiotic, then standardize drug class

    if (antibiotic_col %in% colnames(geno_table_label) | class_col %in% colnames(geno_table_label)) {
      
      # rows where Antibiotic exists
      df_antibiotic <- geno_table_label %>%
        filter(!is.na(!!sym(antibiotic_col)) & !!sym(antibiotic_col) != "") %>%
        separate_longer_delim(!!sym(antibiotic_col), delim = "; ") %>%
        left_join(rgi_drugs, by = setNames("RGI_DrugClassAgent", antibiotic_col)) %>%
        rename(drug_internal = drug_agent) %>%
        rename(drug_class_internal = drug_class) %>%
        mutate(
          drug_to_parse = if_else(!is.na(drug_internal), NA, !!sym(antibiotic_col))
        ) %>%
        mutate(
          drug_agent = AMR::as.ab(drug_to_parse),
          drug_class_agent = AMR::ab_group(drug_to_parse)
        ) %>%
        mutate(
          drug_agent = coalesce(drug_internal, as.character(drug_agent)),
          drug_class = coalesce(drug_class_internal, as.character(drug_class_agent))
        )
      
      # rows where Antibiotic is NA
      df_drugclass <- geno_table_label %>%
        filter(is.na(!!sym(antibiotic_col)) | !!sym(antibiotic_col) == "") %>%
        mutate(
          !!sym(class_col) := if_else(
            (is.na(!!sym(class_col)) | !!sym(class_col) == "") &
              grepl("antibiotic efflux", `Resistance Mechanism`, ignore.case = TRUE),
            "antibiotic efflux",
            !!sym(class_col)
          )
        ) %>%
        filter(!is.na(!!sym(class_col)) & !!sym(class_col) != "") %>%
        separate_longer_delim(!!sym(class_col), delim = "; ") %>%
        left_join(
          rgi_drugs %>% select(RGI_DrugClassAgent, drug_class),
          by = setNames("RGI_DrugClassAgent", class_col)
        )

    # recombine
    geno_table_label_ab <- bind_rows(df_antibiotic, df_drugclass) %>%
      select(-drug_internal, -drug_class_internal, -drug_to_parse, -drug_class_agent) %>%
      dplyr::relocate(any_of(c("id", "marker", "mutation", "drug_agent", "drug_class", "variation type", "marker.label")), .before = dplyr::everything())
    
    } else{
      stop(paste("Input file lacks the expected column: Antibiotic OR `Drug Class` OR `Resistance Mechanism`\n"))
      }

    return(geno_table_label_ab)
}

