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

#' Compare Genotype and Phenotype Data by Sample ID
#'
#' This function compares genotype (`geno_data`) and phenotype (`pheno_data`) datasets based on their sample IDs. It identifies sample IDs that are unique to either dataset, as well as those that overlap. Optionally, it renames the sample ID columns in both datasets for consistency.
#' @param geno_data A data frame containing the genotype data. The first column (or the column specified by `geno_sample_col`) should represent sample IDs.
#' @param pheno_data A data frame containing the phenotype data. The first column (or the column specified by `pheno_sample_col`) should represent sample IDs.
#' @param geno_sample_col A string specifying the name of the column in `geno_data` that contains the sample IDs. Defaults to the first column of `geno_data`.
#' @param pheno_sample_col A string specifying the name of the column in `pheno_data` that contains the sample IDs. Defaults to the first column of `pheno_data`.
#' @param rename_id_cols A logical value indicating whether to rename the sample ID columns in both `geno_data` and `pheno_data` to "id". Defaults to `FALSE`.
#' @importFrom dplyr any_of filter pull rename
#' @return A list containing:
#' - `pheno_unique`: A vector of sample IDs that are unique to the phenotype dataset.
#' - `geno_unique`: A vector of sample IDs that are unique to the genotype dataset.
#' - `overlap_ids`: A vector of sample IDs that are common to both datasets.
#' - `geno_matched`: A data frame of the genotype data filtered to only include the samples with matching IDs.
#' - `pheno_matched`: A data frame of the phenotype data filtered to only include the samples with matching IDs.
#' @export
#' @examples
#' geno_table <- import_amrfp(ecoli_geno_raw)
#' head(ecoli_ast)
#' result <- compare_geno_pheno_id(geno_table,
#'   ecoli_ast,
#'   geno_sample_col = "id",
#'   pheno_sample_col = "id"
#' )
#' print(result$pheno_unique)
#' print(result$geno_unique)
#' print(result$overlap_ids)
#' print(result$geno_matched)
#' print(result$pheno_matched)
compare_geno_pheno_id <- function(geno_data, pheno_data, geno_sample_col = NULL, pheno_sample_col = NULL, rename_id_cols = FALSE) {
  if (is.null(geno_sample_col)) {
    geno_sample_col <- colnames(geno_data)[1]
  }
  if (is.null(pheno_sample_col)) {
    pheno_sample_col <- colnames(pheno_data)[1]
  }

  # get sample ids from both files, extract as vector
  geno_sample_ids <- unique(pull(geno_data, geno_sample_col))
  pheno_sample_ids <- unique(pull(pheno_data, pheno_sample_col))

  # compare the two
  # Unique values in pheno and geno
  unique_in_pheno <- setdiff(pheno_sample_ids, geno_sample_ids)
  unique_in_geno <- setdiff(geno_sample_ids, pheno_sample_ids)

  # Overlapping values
  overlapping_values <- intersect(pheno_sample_ids, geno_sample_ids)

  geno_matched <- geno_data %>% filter(get(geno_sample_col) %in% overlapping_values)
  pheno_matched <- pheno_data %>% filter(get(pheno_sample_col) %in% overlapping_values)

  if (rename_id_cols) {
    geno_matched <- geno_matched %>% rename(id = any_of(geno_sample_col))
    pheno_matched <- pheno_matched %>% rename(id = any_of(pheno_sample_col))
  }

  return(list(
    pheno_unique = unique_in_pheno,
    geno_unique = unique_in_geno,
    overlap_ids = overlapping_values,
    geno_matched = geno_matched,
    pheno_matched = pheno_matched
  ))
}
