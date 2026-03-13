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
#' @param kleborate_class_table A tibble containing a reference table mapping Kleborate drug class column names (`Kleborate_Class`) to standardised drug classes (`drug_class`). Defaults to `kleborate_classes`, which is provided internally.
#' @param hgvs Logical indicating whether to expect mutations in HGVS format (used in Kleborate releases since v3.1.3). Default `TRUE`, which expects mutations formatted as e.g. "GyrA:p.S83F". Set to `FALSE` if your results were generated using older versions where mutations were formatted as e.g. "GyrA_83F".
#' @importFrom dplyr filter left_join mutate select bind_rows rename_with join_by
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
#' # example Kleborate data from EUSCAPE project
#' kleborate_raw
#'
#' # import first few rows of this data frame and parse it to standard genotype table format
#' kleborate_geno <- import_kleborate(kleborate_raw %>% head(n = 10))
#'
#' # parse the output of an older version of Kleborate (v3.1.3) before HGVS 
#' syntax was introduced for mutations
#' kleborate_geno <- import_kleborate(kleborate_raw_v313 %>% head(n = 10), hgvs = FALSE)
import_kleborate <- function(input_table,
                             sample_col = "strain",
                             kleborate_class_table = kleborate_classes,
                             hgvs = TRUE) {
  in_table <- process_input(input_table)

  geno_table <- in_table %>%
    select(any_of(c(sample_col, kleborate_class_table$Kleborate_Class))) %>%
    pivot_longer(-strain, names_to = "Kleborate_Class", values_to = "marker") %>%
    filter(marker != "-") %>%
    tidyr::separate_longer_delim(marker, delim = ";") %>%
    left_join(kleborate_class_table, by = join_by(Kleborate_Class)) %>%
    mutate(marker = str_remove_all(marker, "\\^"))

  # version-specific processing

  if (hgvs) { # newer versions use HGVS nomenclature (e.g. [gene]_:p.[mutation])
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
      relocate(Kleborate_Class, .after = "variation type") %>%
      mutate(drug_agent = NA)
  } else { # older versions use informal nomenclature (e.g. [gene]-[mutation], [gene]-X%, OmpK36GD)
    geno_table <- geno_table %>%
      mutate(`variation type` = case_when(
        grepl("-[0-9]+%$", marker) ~ "Inactivating mutation detected",
        grepl("_mut", Kleborate_Class) & grepl("_[a-z][0-9]+[a-z]$", marker) ~ "Nucleotide variant detected", # Omp mutation c>t
        grepl("_mut", Kleborate_Class) & grepl("-[0-9]+[A-Z]$", marker) ~ "Protein variant detected", # gyrA/parC mutations
        grepl("_mut", Kleborate_Class) & marker %in% c("OmpK36GD", "OmpK36TD") ~ "Protein variant detected",
        TRUE ~ "Gene presence detected"
      )) %>%
      mutate(
        gene = if_else(grepl("%", marker), sub("-.*", "", marker), marker),
        mutation = if_else(grepl("%", marker), sub(".*-", "", marker), NA)
      ) %>%
      mutate(marker.label = if_else(`variation type` == "Inactivating mutation detected",
        paste0(gene, ":-"),
        marker
      )) %>%
      relocate(Kleborate_Class, .after = "variation type") %>%
      mutate(drug_agent = NA)
  }

  return(geno_table)
}
