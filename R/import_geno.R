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

#' Import and Process AMRFinderPlus Results
#'
#' This function imports and processes AMRFinderPlus results, extracting antimicrobial resistance (AMR) elements and mapping them to standardised antibiotic names and drug classes. The function also converts gene symbols to a harmonised format and ensures compatibility with the AMR package.
#' @param input_table A character string specifying a dataframe or path to the AMRFinderPlus results table (TSV format).
#' @param sample_col A character string specifying the column that identifies samples in the dataset (default "`Name`").
#' @param amrfp_drugs A tibble containing a reference table mapping AMRFinderPlus subclasses (`AMRFP_Subclass`) to standardised drug classes (`drug_class`). Defaults to `amrfp_drugs_table`, which is provided internally.
#' @param element_symbol_col Optional character string specifying the column containing gene or element symbols if non-standard column names are used.
#' @param element_type_col Optional character string specifying the column indicating element type (e.g. AMR).
#' @param element_subtype_col Character string specifying the column used to detect mutation subtypes.
#' @param method_col Character string specifying the AMRFinderPlus method column.
#' @param node_col Character string specifying the hierarchy node column.
#' @param subclass_col Character string specifying the AMRFinderPlus subclass column.
#' @param class_col Character string specifying the AMRFinderPlus class column.
#' @importFrom AMR as.ab ab_group
#' @importFrom dplyr all_of everything filter left_join mutate select case_when if_else relocate any_of
#' @importFrom tibble tibble add_column
#' @importFrom tidyr separate_longer_delim separate
#' @importFrom stringr str_match
#' @importFrom rlang sym
#' @importFrom purrr map_chr
#' @return A data frame with the processed genotype data, with harmonised gene names, mapped drug agents, and drug classes which can be used for other functions of the ARMgen package. The output retains the original columns from the AMRFinderPlus table along with the newly mapped variables:
#' - `id`: The sample identifier (`character`).
#' - `marker`: The name of the genotype marker as it appears in the input (e.g. `gyrA_S83F`) (`character`).
#' - `gene`: The gene identifier (`character`).
#' - `mutation`: The mutation detected within the gene, converted to [HGVS nomenclature](https://hgvs-nomenclature.org/stable/) syntax (e.g. `Ser83Phe`) (`character`).
#' - `node`: The node in the NCBI reference gene hierarchy corresponding to the gene (`character`).
#' - `drug_class`: Name of the antibiotic group associated with the genotype marker, compatible with AMR pkg (`character`).
#' - `drug_agent`: Name of the specific antibiotic agent associated with the genotype marker, compatible with AMR pkg (`ab`). Value `NA` is assigned when the markers are annotated with a class only and not a specific antibiotic.
#' - `variation type`: Type of variation, e.g. `Gene presence detected`, `Protein variant detected`, `Nucleotide variant detected`, `Inactivating mutation detected`, `Promoter variant detected`.
#' ... Other fields specific to the input file
#' @details
#' The function performs the following steps:
#' - Reads the AMRFinderPlus output table.
#' - Filters the data to only include AMR elements.
#' - Converts gene symbols to a harmonised format.
#' - Splits multiple subclass annotations into separate rows.
#' - Maps AMRFinderPlus subclasses to standardised drug class names recognised by the AMR pkg.
#' This processing ensures compatibility with downstream AMRgen analysis workflows.
#' @export
#' @examples
#' \dontrun{
#' # small example E. coli AMRFinderPlus data
#' data(ecoli_geno_raw)
#' ecoli_geno_raw
#'
#' # import first few rows of this data frame and parse it as AMRfp data
#' geno <- import_amrfp(ecoli_geno_raw %>% head(n = 10), "Name")
#' geno
#' }
import_amrfp <- function(input_table,
                         sample_col = "Name",
                         element_symbol_col = NULL,
                         element_type_col = NULL,
                         element_subtype_col = "Element subtype",
                         method_col = "Method",
                         node_col = "Hierarchy node",
                         subclass_col = "Subclass",
                         class_col = "Class",
                         amrfp_drugs = amrfp_drugs_table) {
  in_table <- process_input(input_table)

  in_table <- in_table %>% rename(id = !!sym(sample_col))

  if ("Element symbol" %in% colnames(in_table)) {
    in_table <- in_table %>% mutate(marker = `Element symbol`)
    element_symbol_col <- "Element symbol"
  } else if ("Gene symbol" %in% colnames(in_table)) {
    in_table <- in_table %>% mutate(marker = `Gene symbol`)
    element_symbol_col <- "Gene symbol"
  } else if (!is.null(element_symbol_col)) {
    if (element_symbol_col %in% colnames(in_table)) {
      in_table <- in_table %>% mutate(marker = get(element_symbol_col))
    } else {
      stop(paste("Input file lacks the expected column:", element_symbol_col, "\n"))
    }
  } else {
    stop("Input file lacks the expected column: 'Element symbol' (v4.0+) or 'Gene symbol' (pre-v4)\n")
  }

  # filter to only include AMR elements
  if ("Element type" %in% colnames(in_table)) {
    in_table <- in_table %>% filter(`Element type` == "AMR")
  } else if ("Type" %in% colnames(in_table)) {
    in_table <- in_table %>% filter(Type == "AMR")
  } else if (!is.null(element_type_col)) {
    if (element_type_col %in% colnames(in_table)) {
      in_table <- in_table %>% filter(get(element_type_col) == "AMR")
    } else {
      message("Input file lacks the expected column: ", element_type_col, ", assuming all rows report AMR markers.")
    }
  } else {
    message("Input file lacks the expected column: 'Type' (v4.0+) or 'Element type' (pre-v4), assuming all rows report AMR markers.")
  }

  # detect variation type and process mutation
  if (method_col %in% colnames(in_table)) {
    in_table_mutation <- in_table %>%
      mutate(`variation type` = case_when(
        !!sym(method_col) == "INTERNAL_STOP" ~ "Inactivating mutation detected",
        grepl("PARTIAL", !!sym(method_col)) ~ "Inactivating mutation detected",
        !!sym(method_col) == "POINTN" ~ "Nucleotide variant detected",
        !!sym(method_col) %in% c("POINTX", "POINTP") ~ "Protein variant detected",
        !!sym(method_col) %in% c("ALLELEP", "ALLELEX", "BLASTP", "BLASTX", "EXACTP", "EXACTX") ~ "Gene presence detected",
        TRUE ~ NA
      )) %>%
      separate(marker, into = c("gene", "mutation"), sep = "_", remove = FALSE, fill = "right") %>%
      mutate(gene = if_else(startsWith(!!sym(method_col), "POINT"), gene, marker)) %>%
      mutate(mutation = if_else(startsWith(!!sym(method_col), "POINT"), convert_mutation(marker, !!sym(method_col)), mutation))
  } else if (element_subtype_col %in% colnames(in_table)) {
    message("Need method column: ", method_col, " to assign variation type.")
    in_table_mutation <- in_table %>%
      separate(!!sym(element_symbol_col), into = c("gene", "mutation"), sep = "_", remove = FALSE, fill = "right") %>%
      mutate(gene = if_else(startsWith(!!sym(element_subtype_col), "POINT"), gene, marker)) %>%
      mutate(mutation = if_else(startsWith(!!sym(element_subtype_col), "POINT"), purrr::map_chr(marker, convert_mutation, NULL), "-"))
  } else {
    message("Need method column: ", method_col, " or element subtype column: ", element_subtype_col, " columns to parse mutations.")
    in_table_mutation <- in_table %>% mutate(gene = NA, mutation = NA)
  }

  # check for nucleotide variants with negative positions, which indicates promoter variants
  if ("variation type" %in% colnames(in_table_mutation)) {
    in_table_mutation <- in_table_mutation %>%
      mutate(`variation type` = case_when(
        mutation == "-" ~ `variation type`,
        startsWith(mutation, "-") ~ "Promoter variant detected",
        TRUE ~ `variation type`
      ))
  }

  # create AMRrules style label with node:mutation
  if (!(node_col %in% colnames(in_table_mutation))) {
    node_col <- "gene"
  }
  in_table_label <- in_table_mutation %>%
    mutate(node = if_else(is.na(!!sym(node_col)), gene, !!sym(node_col)))

  if (element_subtype_col %in% colnames(in_table_mutation)) {
    in_table_label <- in_table_label %>%
      mutate(marker.label = if_else(!!sym(element_subtype_col) == "POINT",
        if_else(!is.na(mutation), paste0(node, ":", mutation), gsub("_", ":", marker)),
        node
      ))
  } else if (method_col %in% colnames(in_table_mutation)) {
    in_table_label <- in_table_label %>%
      mutate(marker.label = if_else(startsWith(!!sym(method_col), "POINT"),
        if_else(!is.na(mutation), paste0(node, ":", mutation), gsub("_", ":", marker)),
        node
      ))
  } else {
    warning(element_subtype_col, " field not present in input file, guessing mutation markers.")
    in_table_label <- in_table_label %>%
      mutate(marker.label = if_else(mutation != "-", paste0(node, ":", mutation), node))
  }

  if ("variation type" %in% colnames(in_table_label)) {
    in_table_label <- in_table_label %>%
      mutate(marker.label = if_else(`variation type` == "Inactivating mutation detected",
        paste0(node, ":-"),
        marker.label
      ))
  }

  # now split the Subclass column on the "/" to make them one per row, to make adding the ab names easier
  if (subclass_col %in% colnames(in_table_label)) {
    in_table_label <- in_table_label %>% separate_longer_delim(!!sym(subclass_col), "/")
  }

  # make two new columns - drug_class and drug_agent, where we control the vocab for the AMRFinderPlus Subclass column
  # into something that is comparable with the drugs and groups in the AMR package

  # first, identify any subclasses we _know_ aren't in the AMR package, using the internal data
  # join introduces these as new drug_class column
  in_table_ab <- in_table_label %>%
    left_join(amrfp_drugs, by = setNames("AMRFP_Subclass", subclass_col)) %>%
    rename(drug_class_internal = drug_class)

  # then for the columns which are NA, we want to use the Subclass col and convert to ab using AMR pkg
  in_table_ab <- in_table_ab %>%
    mutate(subclass_to_parse = if_else(!is.na(drug_class_internal), NA, !!sym(subclass_col))) %>% # create clean vector of only those subclasses we want to parse with AMR pkg functions
    mutate(drug_agent = AMR::as.ab(subclass_to_parse)) %>%
    mutate(drug_class_from_agent = AMR::ab_group(subclass_to_parse)) %>%
    mutate(drug_class = coalesce(drug_class_internal, drug_class_from_agent)) %>%
    select(-drug_class_from_agent, -drug_class_internal) %>%
    dplyr::relocate(any_of(c("id", "marker", "gene", "mutation", "drug_agent", "drug_class", "variation type", "node", "marker.label")), .before = dplyr::everything())

  return(in_table_ab)
}


#' Import EBI-processed AMRFinderPlus Genotypes from FTP
#'
#' This function imports processed EBI-processed AMRFinderPlus genotyping results. The expected input is genotype data retrieved from the [EBI AMR Portal FTP site](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) either directly or via the function [download_ebi()].
#' Note that files downloaded from the [EBI AMR Portal web browser](https://www.ebi.ac.uk/amr/data/?view=predictions) are formatted differently and can be imported using [import_amrfp_ebi_web].
#'
#' These data are pre-processed by EBI to match NCBI class/subclass to CARD's antibiotic resistance ontology (ARO), however for consistency this function will re-process the data to generate `drug_agent` and `drug_class` fields consistent with the [import_amrfp()] function (the EBI fields `antibiotic*` are also retained).
#' Note several AMRFinderPlus fields are excluded from EBI files, including hierarchy node, method, percent identity and coverage; therefore unlike the [import_amrfp()] function, this function cannot assign `variation type` or `node`.
#' @param input_table R object or file path for the input EBI genotype table (R object, or file path to a TSV or CSV file).
#' @return A data frame with the processed genotype data, with harmonised gene names, mapped drug agents, and drug classes which can be used for other functions of the ARMgen package. The output retains the original columns from the AMRFinderPlus table along with the newly mapped variables:
#' - `id`: The sample identifier (`character`).
#' - `marker`: The name of the genotype marker as it appears in the input (e.g. `gyrA_S83F`) (`character`).
#' - `gene`: The gene identifier (`character`).
#' - `mutation`: The mutation detected within the gene, converted to [HGVS nomenclature](https://hgvs-nomenclature.org/stable/) syntax (e.g. `Ser83Phe`) (`character`).
#' - `drug_class`: Name of the antibiotic group associated with the genotype marker, compatible with AMR pkg (`character`).
#' - `drug_agent`: Name of the specific antibiotic agent associated with the genotype marker, compatible with AMR pkg (`ab`). Value `NA` is assigned when the markers are annotated with a class only and not a specific antibiotic.
#' ... Other fields specific to the input file
#' @details
#' The function performs the following steps:
#' - Reads the EBI-processed genotype table.
#' - Maps AMRFinderPlus subclasses to standardised drug agent and drug class
#'    names using `amrfp_drugs` (EBI-mappings are retained in `antibiotic*` fields.)
#' - Converts drug agent names to the `ab` class from the AMR package.
#' This processing ensures compatibility with downstream AMR analysis workflows.
#' @export
#' @examples
#' \dontrun{
#' # Download quinolone-related genotype data for E. coli, from EBI
#' ebi_geno_raw <- download_ebi(
#'   data = "genotype", species = "Escherichia coli",
#'   geno_subclass = "QUINOLONE"
#' )
#'
#' # Format the file for import
#' ebi_geno <- import_amrfp_ebi_ftp(ebi_geno_raw)
#' }
import_amrfp_ebi_ftp <- function(input_table) {
  input_table <- process_input(input_table) %>%
    rename(id2 = id) # to avoid clash when creating id from BioSample_ID via import_amrfp

  input_table <- import_amrfp(input_table,
    sample_col = "BioSample_ID",
    element_symbol_col = "amr_element_symbol",
    element_type_col = "element_type",
    element_subtype_col = "element_subtype",
    subclass_col = "subclass",
    class_col = "class",
    amrfp_drugs = amrfp_drugs_table
  ) %>%
    select(-any_of("node"))

  return(input_table)
}


#' Import EBI-processed AMRFinderPlus Genotypes from Web
#'
#' This function imports EBI-processed AMRFinderPlus genotyping results. The expected input is genotype data downloaded from the [EBI AMR Portal web browser](https://www.ebi.ac.uk/amr/data/?view=predictions).
#' Note that files downloaded from the [EBI AMR Portal FTP site](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/), either directly or via the function [download_ebi()], are formatted differently and can be imported using [import_amrfp_ebi_ftp].
#'
#' These data are pre-processed by EBI to match NCBI class/subclass to CARD's antibiotic resistance ontology (ARO), however for consistency this function will re-process the data to generate `drug_agent` and `drug_class` fields consistent with the [import_amrfp()] function (the EBI fields `antibiotic*` are also retained).
#' Note several AMRFinderPlus fields are excluded from EBI files, including hierarchy node, method, percent identity and coverage; therefore unlike the [import_amrfp()] function, this function cannot assign `variation type` or `node`.
#' @param input_table R object or file path for the input EBI genotype table (R object, or file path to a TSV or CSV file).
#' @return A data frame with the processed genotype data, with harmonised gene names, mapped drug agents, and drug classes which can be used for other functions of the ARMgen package. The output retains the original columns from the AMRFinderPlus table along with the newly mapped variables:
#' - `id`: The sample identifier (`character`).
#' - `marker`: The name of the genotype marker as it appears in the input (e.g. `gyrA_S83F`) (`character`).
#' - `gene`: The gene identifier (`character`).
#' - `mutation`: The mutation detected within the gene, converted to [HGVS nomenclature](https://hgvs-nomenclature.org/stable/) syntax (e.g. `Ser83Phe`) (`character`).
#' - `drug_class`: Name of the antibiotic group associated with the genotype marker, compatible with AMR pkg (`character`).
#' - `drug_agent`: Name of the specific antibiotic agent associated with the genotype marker, compatible with AMR pkg (`ab`). Value `NA` is assigned when the markers are annotated with a class only and not a specific antibiotic.
#' ... Other fields specific to the input file
#' @details
#' The function performs the following steps:
#' - Reads the EBI-processed genotype table.
#' - Maps AMRFinderPlus subclasses to standardised drug agent and drug class
#'    names using `amrfp_drugs` (EBI-mappings are retained in `antibiotic*` fields.)
#' - Converts drug agent names to the `ab` class from the AMR package.
#' This processing ensures compatibility with downstream AMR analysis workflows.
#' @export
#' @examples
#' \dontrun{
#' # Download data from EBI web portal and import the file
#' ebi_geno <- import_amrfp_ebi_web("amr_records.csv")
#' }
import_amrfp_ebi_web <- function(input_table) {
  input_table <- process_input(input_table)

  colnames(input_table) <- sub("genotype-", "", colnames(input_table))

  input_table <- import_amrfp(input_table,
    sample_col = "BioSample_ID",
    element_symbol_col = "amr_element_symbol",
    element_type_col = "element_type",
    element_subtype_col = "element_subtype",
    subclass_col = "subclass",
    class_col = "class",
    amrfp_drugs = amrfp_drugs_table
  ) %>%
    select(-any_of("node"))

  return(input_table)
}


#' Import EBI-processed AMRFinderPlus Genotypes
#'
#' This function imports EBI-processed AMRFinderPlus genotyping results. The expected input is genotype data downloaded from the [EBI AMR Portal web browser](https://www.ebi.ac.uk/amr/data/?view=predictions), or the [EBI AMR Portal FTP site](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) either directly or via the function [download_ebi()].
#'
#' These data are pre-processed by EBI to match NCBI class/subclass to CARD's antibiotic resistance ontology (ARO), however for consistency this function will re-process the data to generate `drug_agent` and `drug_class` fields consistent with the [import_amrfp()] function (the EBI fields `antibiotic*` are also retained).
#' Note several AMRFinderPlus fields are excluded from EBI files, including hierarchy node, method, percent identity and coverage; therefore unlike the [import_amrfp()] function, this function cannot assign `variation type` or `node`.
#' @param input_table R object or file path for the input EBI genotype table (R object, or file path to a TSV or CSV file).
#' @param web Logical indicating whether the data is from the web portal (default `FALSE`). If `FALSE` input is assumed to be from FTP or [download_ebi].
#' @return A data frame with the processed genotype data, with harmonised gene names, mapped drug agents, and drug classes which can be used for other functions of the ARMgen package. The output retains the original columns from the AMRFinderPlus table along with the newly mapped variables:
#' - `id`: The sample identifier (`character`).
#' - `marker`: The name of the genotype marker as it appears in the input (e.g. `gyrA_S83F`) (`character`).
#' - `gene`: The gene identifier (`character`).
#' - `mutation`: The mutation detected within the gene, converted to [HGVS nomenclature](https://hgvs-nomenclature.org/stable/) syntax (e.g. `Ser83Phe`) (`character`).
#' - `drug_class`: Name of the antibiotic group associated with the genotype marker, compatible with AMR pkg (`character`).
#' - `drug_agent`: Name of the specific antibiotic agent associated with the genotype marker, compatible with AMR pkg (`ab`). Value `NA` is assigned when the markers are annotated with a class only and not a specific antibiotic.
#' ... Other fields specific to the input file
#' @details
#' The function performs the following steps:
#' - Reads the EBI-processed genotype table.
#' - Maps AMRFinderPlus subclasses to standardised drug agent and drug class
#'    names using `amrfp_drugs` (EBI-mappings are retained in `antibiotic*` fields.)
#' - Converts drug agent names to the `ab` class from the AMR package.
#' This processing ensures compatibility with downstream AMR analysis workflows.
#' @export
#' @examples
#' \dontrun{
#' # Download quinolone-related genotype data for E. coli, from EBI
#' ebi_geno_raw <- download_ebi(
#'   data = "genotype", species = "Escherichia coli",
#'   geno_subclass = "QUINOLONE"
#' )
#'
#' # Format the file for import
#' ebi_geno <- import_amrfp_ebi(ebi_geno_raw)
#'
#' # Download data from EBI web portal and import the file
#' ebi_geno_from_web <- import_amrfp_ebi("amr_records.csv", web = TRUE)
#' }
import_amrfp_ebi <- function(input_table, web = FALSE) {
  input_table <- process_input(input_table)

  if (web) {
    import_amrfp_ebi_web(input_table)
  } else {
    import_amrfp_ebi_ftp(input_table)
  }
}

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

#' Import and Process Resistance Gene Identifier (RGI) Results
#'
#' This function imports and processes genotyping results from the Resistance Gene Identifier (RGI, <https://github.com/arpcard/rgi>), extracting antimicrobial resistance determinants and mapping them to standardised drug classes/antibiotics.
#' @param input_table A character string specifying a dataframe or path to the RGI results table (TSV format).
#' @param orf_id_col A character string specifying the column that identifies open reading frame ID (ORF_ID) in the dataset (default `ORF_ID`). This column includes the sample ID and the contig / genomic location and is a default output of RGI.
#' @param sample_id_sep A character string specifying the separator by which the sample ID is separated from the remaining text in `ORF_ID` (Default: `.fasta.txt:`) . For example: in the `ORF_ID` column, "SAMEA3498968.fasta.txt:1_96 # 109511 # 110635....", the sample ID separator is `.fasta.txt:`.
#' @param model_col A character string specifying the column that identifies model type identified by RGI (default `Model_type`).
#' @param antibiotic_col Character string specifying the antibiotic column (default `Antibiotic`).
#' @param class_col Character string specifying the drug class column (default `Drug Class`).
#' @param exclude_loose Logical indicating whether to exclude Loose hits (AMR markers that fall below a curated bitscore cutoff as defined by CARD/RGI). Default `TRUE`, which excludes Loose hits.
#' @param rgi_short_name A tibble containing a reference table mapping model IDs (from CARD/RGI) to shortened model names as provided by CARD (<https://card.mcmaster.ca/download> in aro_index.tsv). Defaults to `rgi_short_name_table`, which is provided internally.
#' @param rgi_drugs A tibble containing a reference table mapping CARD drug class / drug agents to standardised drug classes/names. Defaults to `rgi_drugs_table`, which is provided internally.
#' @param samples_no_amr A vector of sample IDs that have no RGI output because there are no AMR markers identified. For example `c("SampleA", "SampleB")`. (default = `NULL`)
#' @importFrom AMR as.ab ab_group
#' @importFrom dplyr filter left_join mutate select bind_rows anti_join
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
#' # import using sample_id_sep=`_genomic.fna.txt:` and include Loose hits
#' rgi <- import_rgi(rgi_raw, sample_id_sep = "_genomic.fna.txt:", exclude_loose = FALSE)
#'
#' # example RGI data from EuSCAPE project (including only Perfect and Strict hits)
#' rgi_EuSCAPE_raw
#'
#' # import using defaults (sample_id_sep=`.fasta.txt:`, exclude_loose = `TRUE`)
#' import_rgi(rgi_EuSCAPE_raw)
import_rgi <- function(input_table,
                       orf_id_col = "ORF_ID",
                       sample_id_sep = ".fasta.txt:",
                       model_col = "Model_type",
                       antibiotic_col = "Antibiotic",
                       class_col = "Drug Class",
                       exclude_loose = TRUE,
                       rgi_short_name = rgi_short_name_table,
                       rgi_drugs = rgi_drugs_table,
                       samples_no_amr = NULL) {
  in_table <- process_input(input_table)

  in_table <- left_join(in_table, rgi_short_name, by = c("Model_ID" = "Model ID")) %>%
    mutate(id = sub(paste0(sample_id_sep, ".*"), "", .data[[orf_id_col]]))

  in_table[(in_table == "n/a") | (in_table == "")] <- NA

  # Exclude Loose parameter
  if (exclude_loose) { # exclude Loose hits
    geno_table <- in_table %>% filter(`Cut_Off` != "Loose")
  } else { # include all hits
    geno_table <- in_table
  }

  # AMR marker name assignment
  if ("CARD.Short.Name" %in% colnames(geno_table)) { # Use CARD Short Name for marker label
    geno_table <- geno_table %>% mutate(marker = `CARD Short Name`)
  } else if ("Best_Hit_ARO" %in% colnames(geno_table)) { # Use Best_Hit_ARO for marker label (noting that the names could be very long)
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
  } else {
    message("Need method column: Model_type")
  }

  # Check mutation column exists and reformat table to long form - one mutation per row
  if ("SNPs_in_Best_Hit_ARO" %in% colnames(geno_table)) {
    geno_table <- geno_table %>%
      rename(mutation = SNPs_in_Best_Hit_ARO) %>%
      separate_longer_delim(mutation, delim = ", ")
  } else {
    stop(paste("Input file lacks the expected column: SNPs_in_Best_Hit_ARO \n"))
  }

  # create AMRrules style label with gene:mutation
  variant_models <- c(
    "protein variant model",
    "protein overexpression model",
    "rRNA gene variant model"
  )

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

    # Recombine df_antibiotic and df_drugclass
    geno_table_label_ab <- bind_rows(df_antibiotic, df_drugclass) %>%
      select(-drug_internal, -drug_class_internal, -drug_to_parse, -drug_class_agent) %>%
      dplyr::relocate(any_of(c("id", "marker", "mutation", "drug_agent", "drug_class", "variation type", "marker.label")), .before = dplyr::everything())
  } else {
    stop(paste("Input file lacks the expected column: Antibiotic OR `Drug Class` OR `Resistance Mechanism`\n"))
  }

  # All samples in input
  all_samples <- in_table %>%
    select(id) %>%
    distinct()

  # Samples with AMR hits
  samples_with_hits <- geno_table %>%
    filter(!is.na(Cut_Off)) %>%
    select(id) %>%
    distinct()

  # Auto-detected missing samples (sample name in ID column, but remaining rows are NA)
  auto_no_amr <- anti_join(all_samples, samples_with_hits, by = "id")

  # Combine with user input
  if (!is.null(samples_no_amr)) {
    user_no_amr <- tibble(id = samples_no_amr)
    combined_no_amr <- bind_rows(auto_no_amr, user_no_amr) %>% distinct()
  } else {
    combined_no_amr <- auto_no_amr
  }

  # Remove duplicates already in final table
  combined_no_amr <- combined_no_amr %>%
    filter(!id %in% geno_table_label_ab$id)

  # Append and add "No AMR markers" in Cut_Off column just as a marker
  if (nrow(combined_no_amr) > 0 && ncol(geno_table_label_ab) > 0) {
    no_amr_full <- geno_table_label_ab[0, ]
    no_amr_full <- no_amr_full[rep(1, nrow(combined_no_amr)), ]
    no_amr_full$id <- combined_no_amr$id

    if ("Cut_Off" %in% colnames(no_amr_full)) {
      no_amr_full$Cut_Off <- "No AMR markers"
    }

    geno_table_label_ab <- bind_rows(geno_table_label_ab, no_amr_full)
  }
  return(geno_table_label_ab)
}

#' Import and Process ABRicate Results
#'
#' This function imports and processes ABRicate results, extracting antimicrobial resistance (AMR) elements and mapping them to standardised antibiotic names and drug classes. Currently supports results generated using the ResFinder database.
#' @param input_table A character string specifying a dataframe or path to the ABRicate results table.
#' @param sample_col A character string specifying the column that identifies samples in the dataset (default `"FILE"`).
#' @param gene_col A character string specifying the column that identifies gene symbols in the dataset (default `"GENE"`).
#' @param product_col A character string specifying the column that identifies product names in the dataset (default `"PRODUCT"`).
#' @param ab_col A character string specifying the column that identifies which drug/s each detected gene is associated with (default `"RESISTANCE"`).
#' @param db A character string specifying which AMR gene database ABRicate was run with (default `"resfinder"`; `"ncbi"` is also supported).
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
#' - Reads the ABRicate output table via the internal `process_input` function.
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
      warning("'db' parameter ", db, " does not match DATABASE field in input file: ", toString(db_value))
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


#' Convert mutation string based on method
#'
#' This function takes a mutation string (e.g., "gene_REF123ALT") and a
#' mutation method, then extracts and converts parts of the mutation string.
#' Specifically designed for use within `dplyr::mutate()`.
#'
#' @param symbol_col A character vector representing the 'Gene symbol'
#'                        column (the mutation string).
#' @param method_col A character vector representing the 'Method' column.
#' @return A character vector containing the formatted mutation strings
#'         (e.g., "Ala123Trp") or NA if not applicable/match.
#' @export
convert_mutation <- function(symbol_col, method_col) {
  # Ensure inputs are treated as vectors. `mutate` will pass them as such.
  # The regex pattern for splitting the mutation string
  regex_pattern <- "^([^_]+)_([A-Za-z]+)(-?)(\\d+)([A-Za-z]+)$"

  # Apply str_match to the entire gene_symbol_col vector.
  # str_match is already vectorized, so it handles all rows at once.
  result_matrix <- stringr::str_match(symbol_col, regex_pattern)

  # Convert the result matrix to a tibble immediately
  extracted_data <- tibble(
    ref = result_matrix[, 3],
    minus = result_matrix[, 4],
    position = result_matrix[, 5],
    alt = result_matrix[, 6]
  )

  # Perform the conditional conversion based on Method and apply convert_aa_code
  if (!is.null(method_col)) {
    new_ref <- ifelse(method_col %in% c("POINTX", "POINTP"), convert_aa_code(extracted_data$ref), extracted_data$ref)
    new_alt <- ifelse(method_col %in% c("POINTX", "POINTP"), convert_aa_code(extracted_data$alt), extracted_data$alt)
  } else {
    new_ref <- ifelse(extracted_data$minus != "-", convert_aa_code(extracted_data$ref), extracted_data$ref)
    new_alt <- ifelse(extracted_data$minus != "-", convert_aa_code(extracted_data$alt), extracted_data$alt)
  }

  # Construct the final mutation string
  # `paste0` is vectorized. We need to handle NA carefully if new_ref, position, or new_alt are NA.
  if (!is.null(method_col)) {
    new_mutation_string <- ifelse(
      method_col %in% c("POINTX", "POINTP"),
      ifelse(!is.na(new_ref) & !is.na(new_alt) & !is.na(extracted_data$position),
        paste0(new_ref, extracted_data$minus, extracted_data$position, new_alt),
        NA_character_ # Return NA if any of the components are NA
      ),
      paste0(extracted_data$minus, extracted_data$position, new_ref, ">", new_alt)
    )
  } else {
    new_mutation_string <- ifelse(
      extracted_data$minus != "-",
      ifelse(!is.na(new_ref) & !is.na(new_alt) & !is.na(extracted_data$position),
        paste0(new_ref, extracted_data$minus, extracted_data$position, new_alt),
        NA_character_ # Return NA if any of the components are NA
      ),
      paste0(extracted_data$minus, extracted_data$position, new_ref, ">", new_alt)
    )
  }

  return(new_mutation_string)
}


check_mixed_case_grepl <- function(s) {
  # Check if string contains at least one uppercase letter
  has_upper <- grepl("[A-Z]", s)
  # Check if string contains at least one lowercase letter
  has_lower <- grepl("[a-z]", s)

  # Return TRUE only if both conditions are met
  return(has_upper && has_lower)
}


#' Convert single-letter amino acid code(s) to three-letter code(s)
#'
#' This function takes a single-letter amino acid code, a vector of single-letter codes,
#' or a string representing a sequence of single-letter codes. It returns the
#' corresponding three-letter code(s), concatenated directly for sequences.
#'
#' @param input_code A character string (e.g., "A", "MAG") or a vector of
#'                   character strings (e.g., c("A", "G", "C")).
#' @importFrom stringr str_split
#' @return A character string (or vector of strings) with the three-letter
#'         amino acid code(s). For multi-character input strings, a single
#'         concatenated string is returned. Returns NA for individual unmatched codes.
#' @examples
#' # Single character input
#' convert_aa_code("A")
#'
#' # Vector of single characters
#' convert_aa_code(c("M", "A", "G", "Z")) # Z will be NA
#'
#' # Multi-character sequence input
#' convert_aa_code("MAG")
#' convert_aa_code("MAGL")
#'
#' @export
convert_aa_code <- function(input_code) {
  aa_mapping <- c(
    A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", E = "Glu",
    Q = "Gln", G = "Gly", H = "His", I = "Ile", L = "Leu", K = "Lys",
    M = "Met", F = "Phe", P = "Pro", S = "Ser", T = "Thr", W = "Trp",
    Y = "Tyr", V = "Val", `*` = "Ter"
  )
  # `sapply` ensures this works correctly if `input_code` is a vector
  # (which is how `mutate` passes a column to the function).
  sapply(input_code, function(single_input_str) {
    # Handle NA inputs gracefully
    if (is.na(single_input_str)) {
      return(NA_character_)
    }

    # replace any 'STOP' codon with "*"
    single_input_str <- gsub("STOP", "*", single_input_str)

    # if there are any lowercase, assume this is already in 3-letter code and return unchanged
    if (check_mixed_case_grepl(single_input_str)) {
      return(single_input_str)
    }

    # Split the input string into individual characters
    chars <- stringr::str_split(single_input_str, pattern = "")[[1]]

    # Perform the lookup for each character
    converted_chars <- aa_mapping[chars]

    # return NA if any letters can't be converted
    if (any(is.na(converted_chars))) {
      return(NA_character_)
    }

    # Concatenate the 3-letter codes directly (no separator)
    return(paste(unname(converted_chars), collapse = ""))
  }, USE.NAMES = FALSE) # USE.NAMES=FALSE prevents sapply from trying to name the output vector
}


#' Import and process antimicrobial genotype data from common sources
#'
#' This function imports AMR genotyping datasets in formats generated by common bioinformatics tools (AMRFinderPlus, ABRicate, Kleborate, CARD RGI) as well as processed AMRFinderPlus downloadable from EBI.
#' Drug/class annotations given for each genotype marker in the input file are parsed to standard antibiotic names and/or antibiotic groups recognised by the AMR pkg, to facilitate extracting relevant genotype markers for comparison to phenotype data for a specific antibiotic (e.g. using `AMRgen` functions [get_binary_matrix()], [ppv()], [amr_upset()] and [amr_logistic()]).
#' @param input A string representing a dataframe, or a path to an input file, containing the phenotype data in a supported format. These files may be downloaded from public sources such EBI or NCBI, or the files may be generated using common bioinformatics software for AMR genotyping.
#' @param format A string indicating the format of the data: `"amrfp"` (default), `"ebi_web"`, `"ebi_ftp"`, `"kleborate"`, `"abricate"`, `"rgi"`. This determines which importer function the data is passed on to for processing (see below).
#' @param ... Format-specific arguments. See
#' - `"amrfp"` : [import_amrfp()] AMRFinderPlus output
#' - `"ebi_web"` : [import_amrfp_ebi_ftp()] EBI-processed AMRFinderPlus results (from [FTP](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) or [download_ebi()])
#' - `"ebi_ftp"` : [import_amrfp_ebi_web()] EBI-processed AMRFinderPlus results, (from [EBI AMR portal](https://ebi.ac.uk/amr))
#' - `"kleborate"` : [import_kleborate()] Kleborate output
#' - `"abricate"` : [import_abricate()] ABRicate output
#' - `"rgi"` : [import_rgi()] CARD RGI output
#' @return A data frame with the processed genotype data, with harmonised gene names, mapped drug agents, and drug classes which can be used for other functions of the ARMgen package:
#' - `id`: The sample identifier (`character`).
#' - `marker`: The name of the genotype marker as it appears in the input (e.g. `gyrA_S83F`) (`character`).
#' - `gene`: The gene identifier (`character`).
#' - `mutation`: The mutation detected within the gene, converted to [HGVS nomenclature](https://hgvs-nomenclature.org/stable/) syntax (e.g. `Ser83Phe`) (`character`).
#' - `node`: (for AMRFinderPlus input only) The node in the NCBI reference gene hierarchy corresponding to the gene (`character`).
#' - `drug_class`: Name of the antibiotic group associated with the genotype marker, compatible with AMR pkg (`character`).
#' - `drug_agent`: Name of the specific antibiotic agent associated with the genotype marker, compatible with AMR pkg (`ab`). Value `NA` is assigned when the markers are annotated with a class only and not a specific antibiotic.
#' - `variation type`: (for AMRFinderPlus, ABRicate, or Kleborate results) Type of variation, e.g. `Gene presence detected`, `Protein variant detected`, `Nucleotide variant detected`, `Inactivating mutation detected`, `Promoter variant detected`.
#' ... Other fields specific to the input file
#' @export
#' @examples
#' # Import AMRFinderPlus data file
#' data(ecoli_geno_raw)
#' head(ecoli_geno_raw)
#' geno <- import_geno(ecoli_geno_raw %>% head(n = 10), format = "amrfp")
#' head(geno)
#'
#' \dontrun{
#' # Import ABRicate results (that were run using the default db, resfinder)
#' abricate_resfinder <- import_geno("path/to/abricate_resfinder.tsv",
#'   format = "abricate"
#' )
#'
#' # Import ABRicate results that were run using an alternative db (ncbi)
#' abricate_ncbi <- import_geno("path/to/abricate_ncbi.tsv",
#'   format = "abricate",
#'   db = "ncbi"
#' )
#'
#' # Import Kleborate results
#' kleborate_geno <- import_geno(kleborate_raw %>% head(n = 10),
#'   format = "kleborate"
#' )
#'
#' # Import Kleborate results run with an older version without HGVS syntax
#' kleborate_old <- import_geno(kleborate_raw_v313 %>% head(n = 10),
#'   format = "kleborate",
#'   hgvs = FALSE
#' )
#'
#' # Import CARD RGI results with default parameters
#' rgi_geno <- import_geno(rgi_EuSCAPE_raw %>% head(n = 10),
#'   format = "rgi"
#' )
#'
#' # Import CARD RGI results with additional options
#' rgi_geno <- import_geno(rgi_raw %>% head(n = 10),
#'   format = "rgi",
#'   sample_id_sep = "_genomic.fna.txt:", exclude_loose = FALSE
#' )
#'
#' # Download quinolone-related genotype data for E. coli, from EBI
#' ebi_geno_raw <- download_ebi(
#'   data = "genotype", species = "Escherichia coli",
#'   geno_subclass = "QUINOLONE"
#' )
#' # import the downloaded data
#' ebi_geno <- import_geno(ebi_geno_raw,
#'   format = "ebi_ftp"
#' )
#'
#' # Download data from EBI web portal manually, and import the file
#' ebi_geno_from_web <- import_geno("amr_records.csv",
#'   format = "ebi_web"
#' )
#'
#' # Download carbapenem-related genotype data for K. pneumoniae, from NCBI
#' ncbi_geno_raw <- query_ncbi_bq_geno(
#'   taxgroup = "Klebsiella pneumoniae",
#'   geno_subclass = "CARBAPENEM"
#' )
#' # import the downloaded data
#' geno <- import_geno(ncbi_geno_raw,
#'   format = "amrfp",
#'   sample_col = "biosample_acc"
#' )
#' }
import_geno <- function(input,
                        format = "amrfp",
                        ...) {
  fun <- get_geno_importer(format)

  dots <- forward_args(fun, ...)

  rlang::exec(fun,
    input = input,
    !!!dots
  )
}

# registry of geno data import functions, to be dispatched via import_geno()
.geno_importers <- list(
  amrfp = "import_amrfp",
  ebi_ftp = "import_amrfp_ebi_web",
  ebi_web = "import_amrfp_ebi_ftp",
  kleborate = "import_kleborate",
  abricate = "import_abricate",
  rgi = "import_rgi"
)

# function to identify the right genotype importer function to dispatch
# case-insensitive to user-supplied parameter string 'format'
get_geno_importer <- function(format) {
  format <- tolower(format)

  fun_name <- .geno_importers[[format]]

  if (is.null(fun_name)) {
    stop(
      "Unknown format: ", format,
      "\nAvailable formats: ",
      paste(names(.geno_importers), collapse = ", "),
      call. = FALSE
    )
  }

  get(fun_name, envir = asNamespace(utils::packageName()))
}
