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
#' This function imports and processes genotyping results from Kleborate (https://github.com/klebgenomics/Kleborate), extracting antimicrobial resistance determinants and mapping them to standardised drug classes.
#' @param input_table A character string specifying a dataframe or path to the Kleborate results table (TSV format).
#' @param sample_col A character string specifying the column that identifies samples in the dataset (default `strain`).
#' @importFrom dplyr filter left_join mutate select bind_rows rename_with
#' @importFrom tidyr separate_longer_delim separate
#' @importFrom rlang sym
#' @importFrom stringr str_remove_all
#' @return A tibble containing the processed AMR determinants and drug classes that is AMRgen compatible.
#' @details
#' The function performs the following steps:
#' - Reads the Kleborate output table.
#' - Transforms Kleborate output into long form (i.e., one AMR determinant per row).
#' - Maps Kleborate drug classes to standardised drug class names.
#' This processing ensures compatibility with downstream AMRgen analysis workflows.
#' @export
#' @examples
#' \dontrun{
#' # example Kleborate data from EUSCAPE project
#' kleborate_raw
#'
#' # import first few rows of this data frame and parse it as AMRfp data
#' kleborate_geno <- import_kleborate(kleborate_raw %>% head(n = 10), "strain")
#' geno
#' }
import_kleborate <- function(input_table,
                             sample_col = "strain") {
  in_table <- process_input(input_table)

  geno_table <- in_table %>%
    select(any_of(c(sample_col, kleborate_classes$Kleborate_Class))) %>%
    pivot_longer(-strain, names_to = "Kleborate_Class", values_to = "marker") %>%
    filter(marker != "-") %>%
    tidyr::separate_longer_delim(marker, delim = ";") %>%
    left_join(kleborate_classes) %>%
    mutate(marker = str_remove_all(marker, "\\^"))

  geno_table <- geno_table %>%
    mutate(`variation type` = case_when(
      grepl("Ter", marker) ~ "Inactivating mutation detected",
      grepl("del", marker) ~ "Inactivating mutation detected",
      grepl("_mut", Kleborate_Class) & grepl(":p.", marker) ~ "Protein variant detected",
      grepl("_mut", Kleborate_Class) & grepl(":c.", marker) ~ "Nucleotide variant detected",
      TRUE ~ "Gene presence detected"
    )) %>%
    separate(marker, into = c("gene", "mutation"), sep = ":", remove = FALSE, fill = "right") %>%
    mutate(marker.label = if_else(`variation type` == "Inactivating mutation detected",
      paste0(gene, ":-"),
      marker
    )) %>%
    relocate(Kleborate_Class, .after = "variation type")

  return(geno_table)
}
