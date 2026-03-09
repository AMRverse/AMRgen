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
#'
#' Import and Process Abricate Results
#'
#' This function imports and processes Abricate results (e.g., using the ResFinder database), extracting antimicrobial resistance (AMR) elements and mapping them to standardised antibiotic names and drug classes. The function also standardises sample names by replacing dashes with underscores to ensure compatibility with other AMRverse functions.
#' @param input_table A character string specifying a dataframe or path to the Abricate results table.
#' @param sample_col A character string specifying the column that identifies samples in the dataset (default `"Name"`).
#' @param database A dataframe or tibble containing a reference table mapping drug agents (`drug_agent`) to standardised drug classes (`drug_class`). Defaults to `resfinder`, which is provided internally by the package. Currently only the resfinder database is supported. 
#' @importFrom AMR as.ab
#' @importFrom dplyr rename mutate filter select distinct left_join relocate any_of everything
#' @importFrom tidyr separate_longer_delim
#' @importFrom rlang := sym .data
#' @return A tibble containing the processed AMR elements, with harmonised gene names, mapped drug agents, and drug classes which can be used for other functions of the ARMgen package.
#' @details
#' The function performs the following steps:
#' - Reads the Abricate output table via the internal `process_input` function.
#' - Standardises the sample column name and format (replacing hyphens with underscores).
#' - Assigns standardised column names for genes, markers, and variation types (defaulting to "Gene presence detected").
#' - Splits multiple resistance annotations (separated by semicolons) into separate rows.
#' - Maps drug agents to standardised drug class names using the provided `database`.
#' - Converts drug agent names to the `"ab"` class from the AMR package.
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'abricate_raw' is a loaded dataframe of Abricate results
#' abricate_processed <- import_abricate(abricate_raw, sample_col = "Name")
#' head(abricate_processed)
#' }
import_abricate <- function(input_table, 
                            sample_col = "Name", 
                            database = resfinder) {
  
  # Use AMRgen's internal function to handle both file paths and dataframes
  in_table <- process_input(input_table)
  
  # 1. Standardise Sample Column Name
  if ("#FILE" %in% colnames(in_table)) {
    in_table <- in_table %>% dplyr::rename(!!sample_col := `#FILE`)
  } else if ("FILE" %in% colnames(in_table)) {
    in_table <- in_table %>% dplyr::rename(!!sample_col := FILE)
  } else if ("X.FILE" %in% colnames(in_table)) {
    in_table <- in_table %>% dplyr::rename(!!sample_col := X.FILE)
  }

  # 2. Process Core Columns
  in_table_processed <- in_table %>%
    dplyr::mutate(
      marker = GENE, 
      gene = PRODUCT,
      `variation type` = "Gene presence detected",
      mutation = NA_character_,
      node = PRODUCT,
      marker.label = PRODUCT
    )
  
  # 3. Expand Drugs (Handle RESISTANCE column)
  in_table_ab <- in_table_processed %>%
    tidyr::separate_longer_delim(RESISTANCE, delim = ";") %>%
    dplyr::mutate(
      drug_agent = trimws(RESISTANCE) 
    ) %>%
    dplyr::filter(drug_agent != "")
  
  # 4. Join with Database to get 'drug_class'
  # This assumes the database has 'drug_agent' and 'drug_class' columns, which is true of the resfinder mapping file we provide. 
  if (!is.null(database)) {
    # Select only relevant columns to avoid joining unnecessary metadata
    # and use distinct() to prevent row duplication if the reference has duplicates
    mapping_clean <- database %>%
      dplyr::select(dplyr::any_of(c("drug_agent", "drug_class"))) %>%
      dplyr::distinct()
    
    # left_join is performed between the parsed drug_agent from Abricate and the resfinder mapping file to retrieve the corresponding drug_class.
    in_table_ab <- in_table_ab %>%
      dplyr::left_join(mapping_clean, by = "drug_agent")
  } else {
    warning("No database provided; drug_class will be NA.")
    in_table_ab <- in_table_ab %>% 
      dplyr::mutate(drug_class = NA_character_)
  }
  
  # 5. Apply AMR parsing using the AMR package
  in_table_ab <- in_table_ab %>% 
    dplyr::mutate(drug_agent = AMR::as.ab(drug_agent))
  
  # 6. Final Formatting (Relocate to match import_amrfp output)
  in_table_ab <- in_table_ab %>%
    dplyr::relocate(dplyr::any_of(c(
      sample_col, 
      "gene", 
      "mutation", 
      "node", 
      "variation type", 
      "marker", 
      "marker.label",
      "antibiotic",
      "drug_agent", 
      "drug_class"
    )), .before = dplyr::everything())
  
  return(in_table_ab)
}