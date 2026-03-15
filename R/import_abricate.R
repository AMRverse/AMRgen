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
#' @param db A character string specifying which AMR gene database Abricate was run with (default `"resfinder"`; `"ncbi"` is also supported).
#' @importFrom AMR as.ab
#' @importFrom dplyr mutate filter relocate any_of everything rename
#' @importFrom tidyr separate_longer_delim
#' @importFrom rlang := sym .data
#' @return A data frame with the processed genotype data, with harmonised gene names, mapped drug agents, and drug classes which can be used for other functions of the ARMgen package:
#' - `id`: The sample identifier (`character`).
#' - `marker`: The name of the genotype marker, as it appears in the `GENE` column of the input file (`character`).
#' - `gene`: The name of the gene product, as it appears in the `PRODUCT` column of the input file (`character`).
#' - `drug_class`: Name of the antibiotic group associated with the genotype marker, compatible with AMR pkg, parsed from the `RESISTANCE` column of the input file which depends on the database that ABRicate was run with (`character`).
#' - `drug_agent`: Name of the specific antibiotic agent associated with the genotype marker, compatible with AMR pkg, parsed from the `RESISTANCE` column of the input file (`ab`). Value `NA` is assigned when the markers are annotated with a class only and not a specific antibiotic.
#' - `variation type`: Type of variation, i.e. `Gene presence detected`, as ABRicate only detects presence/absence of genes in the query database.
#' ... Other fields specific to the input file
#' @details
#' The function performs the following steps:
#' - Reads the Abricate output table via the internal `process_input` function.
#' - Standardises the sample column name 'id'.
#' - Assigns standardised column names for genes, markers, and sets variation type to "Gene presence detected".
#' - Splits multiple resistance annotations (separated by semicolons) into separate rows.
#' - Converts drug agent names and classes to terms recognised by the AMR package.
#' @export
#' @examples
#' \dontrun{
#' geno_table <- import_abricate("path/to/abricate_resfinder.tsv")
#'
#' geno_table2 <- import_abricate("path/to/abricate_ncbi.tsv", db = "ncbi")
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
  ## TO CHECK: this is coded for resfinder results, where classes are separated by ';'
  ## if using db=ncbi, does this column include multiple subclasses separated by "/" as done in AMRfp?
  in_table <- in_table %>%
    tidyr::separate_longer_delim(!!sym(ab_col), delim = ";") %>%
    dplyr::mutate(!!sym(ab_col) := trimws(!!sym(ab_col))) %>%
    dplyr::filter(!!sym(ab_col) != "")

  if ("DATABASE" %in% colnames(in_table)) {
    db_value <- unique(in_table$DATABASE)
    if (db_value != db) {
      message(paste0("Warning, 'db' parameter ", db, " does not match DATABASE field in input file: ", paste0(db_value, collapse = ", ")))
    }
  }


  # Parse RESISTANCE column values to standard antibiotic and class names used in AMR pkg
  if (db == "ncbi") {
    # first, identify any subclasses we _know_ aren't in the AMR package, using the internal data
    # join introduces these as new drug_class column
    in_table <- in_table %>%
      left_join(amrfp_drugs_table, by = setNames("AMRFP_Subclass", ab_col)) %>%
      rename(drug_class_internal = drug_class)

    # then for the columns which are NA, we want to use the Subclass col and convert to ab using AMR pkg
    in_table <- in_table %>%
      mutate(subclass_to_parse = if_else(!is.na(drug_class_internal), NA, !!sym(ab_col))) %>% # create clean vector of only those subclasses we want to parse with AMR pkg functions
      mutate(drug_agent = AMR::as.ab(subclass_to_parse)) %>%
      mutate(drug_class_from_agent = AMR::ab_group(subclass_to_parse)) %>%
      mutate(drug_class = coalesce(drug_class_internal, drug_class_from_agent))
  } else { # parse drugs directly with AMR package; this works for resfinder
    in_table <- in_table %>%
      mutate(drug_agent = AMR::as.ab(!!sym(ab_col))) %>%
      mutate(drug_class = AMR::ab_group(drug_agent))

    # Harmonise outliers match how we parse NCBI subclass
    in_table <- in_table %>%
      mutate(drug_class = case_when(
        drug_agent %in% c("Sulfamethoxazole", "Sulfathiazole") ~ "Sulfonamides",
        drug_class %in% c("Penicillins", "Aminopenicillins", "Ureidopenicillins", "Monobactams") ~ "Beta-lactams",
        drug_class == "Fluoroquinolones" ~ "Quinolones",
        TRUE ~ drug_class
      ))
  }

  # Move standard AMRgen genotype table cols to the start for visibility
  in_table <- in_table %>%
    dplyr::relocate(dplyr::any_of(c(
      "id",
      "marker",
      "gene",
      "mutation",
      "drug_agent",
      "drug_class",
      "variation type"
    )), .before = dplyr::everything())

  return(in_table)
}
