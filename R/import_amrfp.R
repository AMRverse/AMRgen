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
#' @param sample_col A character string specifying the column that identifies samples in the dataset (default `Name`).
#' @param amrfp_drugs A tibble containing a reference table mapping AMRFinderPlus subclasses (`AFP_Subclass`) to standardised drug agents (`drug_agent`) and drug classes (`drug_class`). Defaults to `amrfp_drugs_table`, which is provided internally.
#' @param element_symbol_col Optional character string specifying the column containing gene or element symbols if non-standard column names are used.
#' @param element_type_col Optional character string specifying the column indicating element type (e.g. AMR).
#' @param element_subtype_col Character string specifying the column used to detect mutation subtypes.
#' @param method_col Character string specifying the AMRFinderPlus method column.
#' @param node_col Character string specifying the hierarchy node column.
#' @param subclass_col Character string specifying the AMRFinderPlus subclass column.
#' @param class_col Character string specifying the AMRFinderPlus class column.
#' @importFrom AMR as.ab
#' @importFrom dplyr all_of everything filter left_join mutate select case_when if_else relocate any_of
#' @importFrom tibble tibble add_column
#' @importFrom tidyr separate_longer_delim separate
#' @importFrom stringr str_match
#' @importFrom rlang sym
#' @importFrom purrr map_chr
#' @return A tibble containing the processed AMR elements, with harmonised gene names, mapped drug agents, and drug classes. The output retains the original columns from the AMRFinderPlus table along with the newly mapped variables.
#' @details
#' The function performs the following steps:
#' - Reads the AMRFinderPlus output table.
#' - Filters the data to only include AMR elements.
#' - Converts gene symbols to a harmonised format.
#' - Splits multiple subclass annotations into separate rows.
#' - Maps AMRFinderPlus subclasses to standardised drug agent and drug class names using `amrfp_drugs`.
#' - Converts drug agent names to the `"ab"` class from the AMR package.
#' This processing ensures compatibility with downstream AMR analysis workflows.
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
      cat(paste0("Input file lacks the expected column: ", element_type_col, ", assuming all rows report AMR markers.\n"))
    }
  } else {
    cat(paste0("Input file lacks the expected column: 'Type' (v4.0+) or 'Element type' (pre-v4), assuming all rows report AMR markers.\n"))
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
    cat("Need method column:", method_col, "to assign variation type.\n")
    in_table_mutation <- in_table %>%
      separate(!!sym(element_symbol_col), into = c("gene", "mutation"), sep = "_", remove = FALSE, fill = "right") %>%
      mutate(gene = if_else(startsWith(!!sym(element_subtype_col), "POINT"), gene, marker)) %>%
      mutate(mutation = if_else(startsWith(!!sym(element_subtype_col), "POINT"), purrr::map_chr(marker, convert_mutation, NULL), "-"))
  } else {
    cat("Need method column:", method_col, "or element subtype column:", element_subtype_col, " columns to parse mutations.\n")
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
    cat(paste("WARNING:", element_subtype_col, "field not present in input file, guessing mutation markers.\n"))
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
  # into something that is comparable with the drugs in the AMR package
  in_table_ab <- in_table_label %>%
    left_join(amrfp_drugs[, c("AFP_Subclass", "drug_agent", "drug_class")], by = setNames("AFP_Subclass", subclass_col)) %>%
    dplyr::relocate(any_of(c(sample_col, "gene", "mutation", "node", "marker", "marker.label", "drug_agent", "drug_class")), .before = dplyr::everything())

  # convert drug_agent into the "ab" class (will leave NAs as is)
  in_table_ab <- in_table_ab %>% mutate(drug_agent = as.ab(drug_agent))

  return(in_table_ab)
}


#' Import EBI-processed AMRFinderPlus Genotypes from FTP
#'
#' This function imports processed EBI-processed AMRFinderPlus genotyping results. The expected input is genotype data retrieved from the [EBI AMR Portal FTP site](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) either directly or via the function [download_ebi()].
#' Note that files downloaded from the [EBI AMR Portal web browser](https://www.ebi.ac.uk/amr/data/?view=predictions) are formatted differently and can be imported using [import_amrfp_ebi_web].
#'
#' These data are pre-processed by EBI to match NCBI class/subclass to CARD's antibiotic resistance ontology (ARO), however for consistency this function will re-process the data to generate `drug_agent` and `drug_class` fields consistent with the [import_amrfp()] function (the EBI fields `antibiotic*` are also retained).
#' Note several AMRfinderplus fields are excluded from EBI files, including hierarchy node, method, percent identity and coverage; therefore unlike the [import_amrfp()] function, this function cannot assign `variation type` or `node`.
#' @param input_table R object or file path for the input EBI genotype table (R object, or file path to a TSV or CSV file).
#' @return A tibble containing the processed AMR elements, with harmonised gene names, drug agents, and drug classes. The output retains the original columns from the AMRFinderPlus table along with the newly mapped variables.
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
  input_table <- process_input(input_table)

  input_table <- import_amrfp(input_table,
    sample_col = "BioSample_ID",
    element_symbol_col = "amr_element_symbol",
    element_type_col = "element_type",
    element_subtype_col = "element_subtype",
    subclass_col = "subclass",
    class_col = "class",
    amrfp_drugs = amrfp_drugs_table
  )

  return(input_table)
}


#' Import EBI-processed AMRFinderPlus Genotypes from Web
#'
#' This function imports EBI-processed AMRFinderPlus genotyping results. The expected input is genotype data downloaded from the [EBI AMR Portal web browser](https://www.ebi.ac.uk/amr/data/?view=predictions).
#' Note that files downloaded from the [EBI AMR Portal FTP site](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/), either directly or via the function [download_ebi()], are formatted differently and can be imported using [import_amrfp_ebi_ftp].
#'
#' These data are pre-processed by EBI to match NCBI class/subclass to CARD's antibiotic resistance ontology (ARO), however for consistency this function will re-process the data to generate `drug_agent` and `drug_class` fields consistent with the [import_amrfp()] function (the EBI fields `antibiotic*` are also retained).
#' Note several AMRfinderplus fields are excluded from EBI files, including hierarchy node, method, percent identity and coverage; therefore unlike the [import_amrfp()] function, this function cannot assign `variation type` or `node`.
#' @param input_table R object or file path for the input EBI genotype table (R object, or file path to a TSV or CSV file).
#' @return A tibble containing the processed AMR elements, with harmonised gene names, drug agents, and drug classes. The output retains the original columns from the AMRFinderPlus table along with the newly mapped variables.
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
  )


  return(input_table)
}


#' Import EBI-processed AMRFinderPlus Genotypes
#'
#' This function imports EBI-processed AMRFinderPlus genotyping results. The expected input is genotype data downloaded from the [EBI AMR Portal web browser](https://www.ebi.ac.uk/amr/data/?view=predictions), or the [EBI AMR Portal FTP site](https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/) either directly or via the function [download_ebi()].
#'
#' These data are pre-processed by EBI to match NCBI class/subclass to CARD's antibiotic resistance ontology (ARO), however for consistency this function will re-process the data to generate `drug_agent` and `drug_class` fields consistent with the [import_amrfp()] function (the EBI fields `antibiotic*` are also retained).
#' Note several AMRfinderplus fields are excluded from EBI files, including hierarchy node, method, percent identity and coverage; therefore unlike the [import_amrfp()] function, this function cannot assign `variation type` or `node`.
#' @param input_table R object or file path for the input EBI genotype table (R object, or file path to a TSV or CSV file).
#' @param web Logical indicating whether the data is from the web portal (default `FALSE`). If `FALSE` input is assumed to be from FTP or [download_ebi].
#' @return A tibble containing the processed AMR elements, with harmonised gene names, drug agents, and drug classes. The output retains the original columns from the AMRFinderPlus table along with the newly mapped variables.
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
