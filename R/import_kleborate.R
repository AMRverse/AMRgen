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
#' @param hgvs Logical indicating whether to expect mutations in [HGVS nomenclature](https://hgvs-nomenclature.org/stable/) syntax (used in Kleborate releases since v3.1.3). Default `TRUE`, which expects mutations formatted as e.g. "GyrA:p.S83F". Set to `FALSE` if your results were generated using older versions where mutations were formatted as e.g. "GyrA_83F".
#' @importFrom dplyr filter left_join mutate select bind_rows rename_with join_by
#' @importFrom tidyr separate_longer_delim separate
#' @importFrom rlang sym
#' @importFrom stringr str_remove_all
#' @return A data frame with the processed genotype data, with harmonised gene names, mapped drug agents, and drug classes which can be used for other functions of the ARMgen package:
#' - `id`: The sample identifier (`character`).
#' - `marker`: The name of the genotype marker as it appears in the input (e.g. `GyrA:p.S83F` for recent versions of Kleborate, or `GyrA-83F` for earlier versions not using [HGVS nomenclature](https://hgvs-nomenclature.org/stable/)) (`character`).
#' - `gene`: The gene identifier (`character`).
#' - `mutation`: The mutation detected within the gene, converted to [HGVS nomenclature](https://hgvs-nomenclature.org/stable/) syntax (e.g. `Ser83Phe`) (`character`).
#' - `drug_class`: Name of the antibiotic group associated with the genotype marker, compatible with AMR pkg, drawn from the Kleborate column in which the marker was reported (`character`).
#' - `drug_agent`: Values are recorded as `NA` as Kleborate doesn't report markers assigned to individual drug level.
#' - `variation type`: Type of variation, e.g. `Gene presence detected`, `Protein variant detected`, `Nucleotide variant detected`, `Inactivating mutation detected`.
#' ... Other fields specific to the input file
#' @details
#' The function performs the following steps:
#' - Reads the Kleborate output table.
#' - Transforms Kleborate output into long form (i.e., one AMR determinant per row).
#' - Maps Kleborate drug classes to standardised drug class names recognised by the AMR pkg.
#' This processing ensures compatibility with downstream AMRgen analysis workflows.
#' @examples
#' # example Kleborate data from EUSCAPE project
#' kleborate_raw
#'
#' # import first few rows of this data frame and parse it to standard genotype table format
#' kleborate_geno <- import_kleborate(kleborate_raw %>% head(n = 10))
#'
#' # parse the output of an older version of Kleborate (v3.1.3) before
#' # HGVS syntax was introduced for mutations
#' kleborate_geno <- import_kleborate(kleborate_raw_v313 %>% head(n = 10), hgvs = FALSE)
#' @export
import_kleborate <- function(input_table,
                             sample_col = "strain",
                             kleborate_class_table = kleborate_classes,
                             hgvs = TRUE) {
  in_table <- process_input(input_table)

  in_table <- in_table %>% rename(id = !!sym(sample_col))

  geno_table <- in_table %>%
    select(any_of(c("id", kleborate_class_table$Kleborate_Class))) %>%
    pivot_longer(-id, names_to = "Kleborate_Class", values_to = "marker") %>%
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
