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
#' This function imports and processes Abricate results, extracting antimicrobial resistance (AMR) elements and mapping them to standardised antibiotic names and drug classes. Currently supports results generated using the ResFinder database.
#' @param input_table A character string specifying a dataframe or path to the Abricate results table.
#' @param sample_col A character string specifying the column that identifies samples in the dataset (default `"FILE"`).
#' @param gene_col A character string specifying the column that identifies gene symbols in the dataset (default `"GENE"`).
#' @param product_col A character string specifying the column that identifies product names in the dataset (default `"PRODUCT"`).
#' @param ab_col A character string specifying the column that identifies which drug/s each detected gene is associated with (default `"RESISTANCE"`).
#' @param db A character string specifying which AMR gene database Abricate was run with (currently only `"resfinder"` is supported).
#' @importFrom AMR as.ab
#' @importFrom dplyr mutate filter relocate any_of everything rename
#' @importFrom tidyr separate_longer_delim
#' @importFrom rlang := sym .data
#' @return A tibble containing the processed AMR elements, with harmonised gene names, mapped drug agents, and drug classes which can be used for other functions of the ARMgen package.
#' @details
#' The function performs the following steps:
#' - Reads the Abricate output table via the internal `process_input` function.
#' - Standardises the sample column name 'id'.
#' - Assigns standardised column names for genes, markers, and sets variation type to "Gene presence detected".
#' - Splits multiple resistance annotations (separated by semicolons) into separate rows.
#' - Converts drug agent names to the `"ab"` class from the AMR package and maps these to classes compatible with the output of [import_amrfp()] and [import_kleborate()].
#' @export
#' @examples
#' \dontrun{
#' geno_table <- import_abricate("path/to/abricate_resfinder.tsv")
#' }
import_abricate <- function(input_table,
                            sample_col = "FILE",
                            gene_col = "GENE",
                            product_col = "PRODUCT",
                            ab_col = "RESISTANCE",
                            db = "resfinder") {
  # note this also strips the # from the start of the file
  in_table <- process_input(input_table)

  in_table <- in_table %>% rename(id = !!sym(sample_col))

  # Process Core Columns
  in_table <- in_table %>%
    dplyr::mutate(
      marker := !!sym(gene_col),
      gene := !!sym(product_col),
      `variation type` = "Gene presence detected",
      mutation = NA_character_
    )

  # Expand Drugs (Handle RESISTANCE column)
  in_table <- in_table %>%
    tidyr::separate_longer_delim(!!sym(ab_col), delim = ";") %>%
    dplyr::mutate(!!sym(ab_col) := trimws(!!sym(ab_col))) %>%
    dplyr::filter(!!sym(ab_col) != "")

  # ResFinder drug names can be parsed directly with AMR package
  # to support output run with other dbs we will need to update this to run through some options
  in_table <- in_table %>%
    mutate(drug_agent = AMR::as.ab(!!sym(ab_col))) %>%
    mutate(drug_class = AMR::ab_group(drug_agent))

  # Harmonise outliers to harmonise classes with how we parse NCBI AMRfinderplus subclass for consistency
  in_table <- in_table %>%
    mutate(drug_class = case_when(
      drug_agent %in% c("Sulfamethoxazole", "Sulfathiazole") ~ "Sulfonamides",
      drug_class %in% c("Penicillins", "Aminopenicillins", "Ureidopenicillins", "Monobactams") ~ "Beta-lactams",
      drug_class == "Fluoroquinolones" ~ "Quinolones",
      TRUE ~ drug_class
    ))

  # Move standard AMRgen genotype table cols to the start for visibility
  in_table <- in_table %>%
    dplyr::relocate(dplyr::any_of(c(
      "id",
      "gene",
      "mutation",
      "variation type",
      "marker",
      "drug_agent",
      "drug_class"
    )), .before = dplyr::everything())

  return(in_table)
}
