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
#'   \item{per_type}{A tibble of unique counts of samples, markers, genes, drugs, and classes per variation type.}
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
                           sample_col = "id",
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
  uniques_pervartype <- NULL
  if (!is.null(variation_col)) {
    if (variation_col %in% colnames(geno_table)) {
      uniques_pervartype <- geno_table %>%
        select(any_of(c(sample_col, marker_col, drug_col, class_col, gene_col, variation_col))) %>%
        group_by(!!sym(variation_col)) %>%
        summarise(across(setdiff(names(.), variation_col), n_distinct))
    }
  }

  # drugs
  drugs <- NULL
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

  markers <- NULL
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
    per_type = uniques_pervartype,
    drugs = drugs,
    markers = markers
  ))
}


#' Summarise a Phenotype Table
#'
#' `summarise_pheno()` computes summary information for a phenotype table.
#'
#' @param pheno_table A tibble or data frame containing phenotype data, in the format output by [import_ast].
#' @param sample_col Character. Name of the column containing sample identifiers. Default is `"id"`.
#' @param drug_col Character. Name of the column containing drug agent identifiers. Default is `"drug_agent"`. If this is of class 'ab' the entries will be annotated with their full antibiotic names, converted using [as.ab]. If this is desired behaviour but the class is not 'ab', set `force_ab=TRUE`.
#' @param mic_col Character. Name of the column containing MIC data. Default is `"mic"`.
#' @param disk_col Character. Name of the column containing drug class identifiers. Default is `"disk"`.
#' @param spp_col Character. Name of the column containing species names Default is `"spp_pheno"`.
#' @param pheno_cols Vector. Vector giving names of columns containing categorical phenotype calls (S/I/R or NWT/WT). Default is any columns beginning with `"pheno"` or `"ecoff"`.
#' @param method_cols Vector. Vector giving names of columns containing method or source information by which to summarise MIC/disk data. Default is `c("method", "platform", "guideline", "source")`.
#' @param force_ab Logical. If `TRUE`, attempts to convert entries in `drug_col` to antibiotic names using [as.ab] even if this column is not of class `"ab"` Default is `FALSE`.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{uniques}{A tibble of the number of unique samples, drugs, organisms, and methods detected in `pheno_table`.}
#'   \item{drugs}{A tibble listing the drugs included in the table, and the associated number of samples with MIC measures, disk measures, neither or both, for each drug and species.}
#'   \item{details}{A tibble listing more details of the methods of assay measurements, per drug and species.}
#'   \item{pheno_counts_list}{A list of tibbles, each corresponding to a unique categorical phenotype column in the input, indicating the counts of each phenotypic category per drug and species.}
#' }
#'
#' @details
#' The function automatically adapts to the presence or absence of columns in `pheno_table`.
#' The `force_ab` parameter allows the addition of full antibiotic names using the `ab_name()` function even when the first column is not recognized as an `"ab"` object.
#' @importFrom dplyr summarise across select count mutate n_distinct
#' @importFrom tidyr pivot_longer pivot_wider
#' @examples
#'
#' summarise_pheno(staph_ast_ebi)
#'
#' summarise_pheno(staph_ast_ebi, pheno_cols = c("pheno_provided", "pheno_clsi", "ecoff"))
#'
#' @export
summarise_pheno <- function(pheno_table,
                            sample_col = "id",
                            drug_col = "drug_agent",
                            mic_col = "mic",
                            disk_col = "disk",
                            spp_col = "spp_pheno",
                            pheno_cols = NULL,
                            method_cols = c("method", "platform", "guideline", "source"),
                            force_ab = FALSE) {
  # uniques per column
  uniques <- pheno_table %>%
    summarise(across(
      any_of(c(sample_col, drug_col, spp_col, method_cols)),
      ~ n_distinct(.x, na.rm = FALSE)
    )) %>%
    tidyr::pivot_longer(everything(),
      names_to = "column",
      values_to = "n_unique"
    )

  # ensure we have columns for mic and disk, set to NA if they were not in the input
  if (is.null(mic_col)) {
    pheno_table <- mutate(mic = NA)
    mic_col <- "mic"
    cat("No MIC data colummn provided\n")
  }
  if (is.null(disk_col)) {
    pheno_table <- mutate(disk = NA)
    disk_col <- "disk"
    cat("No disk data colummn provided\n")
  }

  pheno_table <- pheno_table %>%
    mutate(INTERNAL_measures = case_when(
      is.na(mic) & is.na(disk) ~ "none",
      is.na(mic) & !is.na(disk) ~ "disk",
      !is.na(mic) & is.na(disk) ~ "mic",
      !is.na(mic) & !is.na(disk) ~ "both"
    ))

  drugs <- pheno_table %>%
    select(drug_col, spp_col, "INTERNAL_measures") %>%
    count(across(everything())) %>%
    tidyr::pivot_wider(names_from = "INTERNAL_measures", values_from = n)

  if (inherits(drugs[[1]], "ab") | force_ab) {
    drugs <- drugs %>%
      mutate(antibiotic_name = ab_name(!!sym(drug_col)), .after = !!sym(drug_col))
  }

  if (spp_col %in% colnames(drugs)) {
    if (inherits(drugs[[spp_col]], "mo")) {
      drugs <- drugs %>%
        mutate(!!sym(spp_col) := mo_name(!!sym(spp_col)))
    }
  }


  # add pheno counts where available
  pheno_counts_list <- list()
  if (is.null(pheno_cols)) {
    pheno_cols_list <- pheno_table %>%
      select(starts_with("pheno"), starts_with("ecoff")) %>%
      colnames()
    cat("No phenotype column names provided via 'pheno_cols'\n")
    cat("These are needed to summarise counts of phenotype category calls per drug.\n")
    cat(paste0("Relevant columns detected in your input table are: c('", paste(pheno_cols_list, collapse = "','"), "')\n"))
  }
  if (!is.null(pheno_cols)) {
    for (pheno_col in pheno_cols) {
      pheno_counts <- pheno_table %>%
        select(any_of(c(drug_col, spp_col, pheno_col))) %>%
        count(across(everything())) %>%
        tidyr::pivot_wider(names_from = pheno_col, values_from = n)
      if (inherits(pheno_counts[[1]], "ab") | force_ab) {
        pheno_counts <- pheno_counts %>%
          mutate(antibiotic_name = ab_name(!!sym(drug_col)), .after = !!sym(drug_col))
      }
      if (spp_col %in% colnames(pheno_counts)) {
        if (inherits(pheno_counts[[spp_col]], "mo")) {
          pheno_counts <- pheno_counts %>%
            mutate(!!sym(spp_col) := mo_name(!!sym(spp_col)))
        }
      }
      pheno_counts_list[[pheno_col]] <- pheno_counts
    }
  }

  if (!is.null(method_cols)) {
    details <- pheno_table %>%
      select(any_of(c(drug_col, spp_col, method_cols, "INTERNAL_measures"))) %>%
      count(across(everything())) %>%
      tidyr::pivot_wider(names_from = "INTERNAL_measures", values_from = n)

    if (inherits(drugs[[1]], "ab") | force_ab) {
      details <- details %>%
        mutate(antibiotic_name = ab_name(!!sym(drug_col)), .after = !!sym(drug_col))
    }

    if (spp_col %in% colnames(details)) {
      if (inherits(details[[spp_col]], "mo")) {
        details <- details %>%
          mutate(!!sym(spp_col) := mo_name(!!sym(spp_col)))
      }
    }
  }

  return(list(
    uniques = uniques,
    drugs = drugs,
    details = details,
    pheno_counts_list = pheno_counts_list
  ))
}

# summarise_geno_pheno(geno_table, pheno_table)
