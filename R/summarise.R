#' Summarise a Genotype Table
#'
#' `summarise_geno()` computes summary information for a genotype table.
#'
#' @param geno_table A tibble or data frame containing genotype data, in the format output by [import_amrfp].
#' @param sample_col Character. Name of the column containing sample identifiers. Default is `"Name"`.
#' @param marker_col Character. Name of the column containing marker identifiers. Default is `"marker"`.
#' @param drug_col Character. Name of the column containing drug agent identifiers. Default is `"drug_agent"`. If this is of class 'ab' the entries will be annotated with their full antibiotic names, converted using [as.ab]. If this is desired behaviour but the class is not 'ab', set `force_ab=TRUE`.
#' @param class_col Character. Name of the column containing drug class identifiers. Default is `"drug_class"`.
#' @param gene_col Character. Name of the column containing gene identifiers. Default is `"gene"`.
#' @param variation_col Character. Name of the column containing variation type identifiers. Default is `"variation type"`.
#' @param force_ab Logical. If `TRUE`, attempts to convert entries in `drug_col` to antibiotic names using [as.ab] even if this column is not of class `"ab"` Default is `FALSE`.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{uniques}{A tibble of the number of unique samples, markers, genes, drugs, classes and variation types detected in `geno_table`.}
#'   \item{pertype}{A tibble of unique counts of samples, markers, genes, drugs, and classes per variation type.}
#'   \item{drugs}{A tibble listing the drugs and/or drug classes represented in the table, and the associated number of unique markers, unique samples, and total hits for each drug/class.}
#'   \item{markers}{A tibble listing the markers represented in the table, and the associated drugs/classes and variation types (if present). Number indicates the count of hits detected per marker.}
#' }
#'
#' @details
#' The function automatically adapts to the presence or absence of columns in `geno_table`.
#' The `force_ab` parameter allows the addition of full antibiotic names using the `ab_name()` function even when the first column is not recognized as an `"ab"` object.
#' @importFrom dplyr summarise across select group_by count rename full_join mutate n_distinct
#' @importFrom tidyr pivot_longer
#' @examples
#'
#' geno_table <- import_amrfp(ecoli_geno_raw)
#' summarise_geno(geno_table)
#'
#' @export
summarise_geno <- function(geno_table,
                           sample_col = "Name",
                           marker_col = "marker",
                           drug_col = "drug_agent",
                           class_col = "drug_class",
                           gene_col = "gene",
                           variation_col = "variation type",
                           force_ab = FALSE) {
  # uniques per column
  uniques <- geno_table %>%
    summarise(across(
      any_of(c(sample_col, marker_col, drug_col, class_col, gene_col, variation_col)),
      ~ n_distinct(.x, na.rm = FALSE)
    )) %>%
    tidyr::pivot_longer(everything(),
      names_to = "column",
      values_to = "n_unique"
    )

  # uniques per variation type
  uniques_pervartype <- geno_table %>%
    select(any_of(c(sample_col, marker_col, drug_col, class_col, gene_col, variation_col))) %>%
    group_by(!!sym(variation_col)) %>%
    summarise(across(setdiff(names(.), variation_col), n_distinct))

  # drugs
  if (drug_col %in% colnames(geno_table)) {
    if (class_col %in% colnames(geno_table)) {
      # agents and classes
      drugs <- geno_table %>%
        count(!!sym(drug_col), !!sym(class_col)) %>%
        rename(hits = n)
      if (sample_col %in% colnames(geno_table)) {
        drugs <- geno_table %>%
          count(!!sym(drug_col), !!sym(class_col), !!sym(sample_col)) %>%
          count(!!sym(drug_col), !!sym(class_col)) %>%
          rename(samples = n) %>%
          full_join(drugs, by = c(drug_col, class_col))
      }
      if (marker_col %in% colnames(geno_table)) {
        drugs <- geno_table %>%
          count(!!sym(drug_col), !!sym(class_col), !!sym(marker_col)) %>%
          count(!!sym(drug_col), !!sym(class_col)) %>%
          rename(markers = n) %>%
          full_join(drugs, by = c(drug_col, class_col))
      }
    } else {
      # only have agents, no classes
      drugs <- geno_table %>%
        count(!!sym(drug_col))
      if (sample_col %in% colnames(geno_table)) {
        drugs <- geno_table %>%
          count(!!sym(drug_col), !!sym(sample_col)) %>%
          count(!!sym(drug_col)) %>%
          rename(samples = n) %>%
          full_join(drugs, by = c(drug_col))
      }
      if (marker_col %in% colnames(geno_table)) {
        drugs <- geno_table %>%
          count(!!sym(drug_col), !!sym(marker_col)) %>%
          count(!!sym(drug_col)) %>%
          rename(markers = n) %>%
          full_join(drugs, by = c(drug_col))
      }
    }
    # add full drug name
    if (inherits(drugs[[1]], "ab") | force_ab) {
      drugs <- drugs %>% mutate(antibiotic = ab_name(!!sym(drug_col)), .after = !!sym(drug_col))
    }
  } else if (class_col %in% colnames(geno_table)) {
    # only have classes
    drugs <- geno_table %>%
      count(!!sym(class_col))
    if (sample_col %in% colnames(geno_table)) {
      drugs <- geno_table %>%
        count(!!sym(class_col), !!sym(sample_col)) %>%
        count(!!sym(class_col)) %>%
        rename(samples = n) %>%
        full_join(drugs, by = c(class_col))
    }
    if (marker_col %in% colnames(geno_table)) {
      drugs <- geno_table %>%
        count(!!sym(class_col), !!sym(marker_col)) %>%
        count(!!sym(class_col)) %>%
        rename(markers = n) %>%
        full_join(drugs, by = c(class_col))
    }
  }

  if (marker_col %in% colnames(geno_table)) {
    markers <- geno_table %>%
      select(any_of(c(marker_col, drug_col, class_col, variation_col))) %>%
      count(across(everything()))
    # add full drug name
    if ((inherits(drugs[[1]], "ab") | force_ab) & (drug_col %in% colnames(markers))) {
      markers <- markers %>%
        mutate(antibiotic = ab_name(!!sym(drug_col)), .after = !!sym(drug_col))
    }
  }

  return(list(
    uniques = uniques,
    pertype = uniques_pervartype,
    drugs = drugs,
    markers = markers
  ))
}

# summarise_pheno(pheno_table)

# summarise_geno_pheno(geno_table, pheno_table)
