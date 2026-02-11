#' Get Clinical Breakpoints for an Antibiotic
#'
#' This function retrieves the clinical breakpoints for a given species, antibiotic, and guideline, from the AMR package.
#' It attempts to find the breakpoints at various taxonomic levels (species, genus, family, order) if no direct match is found.
#'
#' @param species A character string representing the species of interest (e.g., "Escherichia coli").
#' @param guide A character string indicating the guideline for breakpoints (default, "EUCAST 2024").
#' @param antibiotic A character string indicating the antibiotic for which to retrieve breakpoints (e.g., "Ciprofloxacin").
#' @param type_filter A character string indicating the type of breakpoints to retrieve (e.g., "human"). Default is "human" which returns human clinical breakpoints, change to "ECOFF" to get the epidemiological cutoff.
#'
#' @return A data frame containing the clinical breakpoints for the specified species, antibiotic, and guideline.
#' If no exact match is found, the function attempts to retrieve breakpoints at the genus, family, and order levels.
#'
#' @examples
#' getBreakpoints("Escherichia coli", "EUCAST 2024", "Ciprofloxacin")
#'
#' @importFrom AMR as.mo as.ab
#' @importFrom dplyr filter
#' @export
getBreakpoints <- function(species, guide = "EUCAST 2024", antibiotic, type_filter = "human") {
  bp <- AMR::clinical_breakpoints %>%
    filter(guideline == guide & mo == AMR::as.mo(species) & ab == AMR::as.ab(antibiotic)) %>%
    filter(type == type_filter)
  if (nrow(bp) == 0) {
    sp_mo <- AMR::as.mo(species)
    this_mo <- AMR::microorganisms %>% filter(mo == sp_mo)
    # try genus
    bp <- AMR::clinical_breakpoints %>%
      filter(guideline == guide & mo == AMR::as.mo(this_mo$genus) & ab == AMR::as.ab(antibiotic)) %>%
      filter(type == type_filter)
    if (nrow(bp) == 0) {
      # try family
      bp <- AMR::clinical_breakpoints %>%
        filter(guideline == guide & mo == AMR::as.mo(this_mo$family) & ab == AMR::as.ab(antibiotic)) %>%
        filter(type == type_filter)
      if (nrow(bp) == 0) {
        # try order
        bp <- AMR::clinical_breakpoints %>%
          filter(guideline == guide & mo == AMR::as.mo(this_mo$order) & ab == AMR::as.ab(antibiotic)) %>%
          filter(type == type_filter)
      }
    }
  }
  return(bp)
}

#' Check and Retrieve Breakpoints for an Antibiotic
#'
#' This function checks the clinical breakpoints for a specified antibiotic and species using the `getBreakpoints` function.
#' It handles cases where multiple breakpoint sites exist and uses the specified site or the one with the highest susceptibility
#' breakpoint.
#'
#' @param species A character string representing the species of interest (e.g., "Escherichia coli").
#' @param guide A character string indicating the guideline for breakpoints (default, "EUCAST 2024").
#' @param antibiotic A character string indicating the antibiotic for which to check breakpoints (e.g., "Ciprofloxacin").
#' @param bp_site A character string specifying the breakpoint site to use (optional). If provided, the function uses this site; otherwise, if different breakpoints are specified for different sites it selects the one with the highest susceptibility breakpoint.
#' @param assay A character string specifying the assay type (either "MIC" or "Disk"). Default is "MIC".
#'
#' @return A list containing:
#' - \code{breakpoint_S}: The susceptibility breakpoint (e.g., MIC or disk size).
#' - \code{breakpoint_R}: The resistance breakpoint (e.g., MIC or disk size).
#' - \code{bp_standard}: The breakpoint site used, if multiple breakpoints are set in the guidelines.
#'
#' @details
#' The function first attempts to retrieve breakpoints using the `getBreakpoints` function. If multiple breakpoint sites are found, it handles the situation by:
#' - Using the specified site if it exists.
#' - Selecting the breakpoint with the highest susceptibility value if the specified site is not found.
#' - Returning a message about the selected site and breakpoint values.
#'
#' @examples
#' checkBreakpoints(
#'   species = "Escherichia coli", guide = "EUCAST 2024",
#'   antibiotic = "Ciprofloxacin", assay = "MIC"
#' )
#'
#' @importFrom dplyr filter arrange first
#' @export
checkBreakpoints <- function(species, guide = "EUCAST 2024", antibiotic, assay = "MIC", bp_site = NULL) {
  breakpoints <- getBreakpoints(species, guide, antibiotic) %>% filter(method == assay)
  if (nrow(breakpoints) == 0) {
    stop(paste("Could not determine", assay, "breakpoints using AMR package, please provide your own breakpoints"))
  } else {
    breakpoint_sites <- unique(breakpoints$site)
    breakpoint_message_multibp <- NA
    bp_standard <- "-"
    # handle multiple breakpoints (e.g. for different conditions)
    if (length(breakpoint_sites) > 1) {
      breakpoint_message_multibp <- paste("NOTE: Multiple breakpoint entries, for different sites:", paste(breakpoint_sites, collapse = "; "))
      if (length(unique(breakpoints$breakpoint_R)) == 1 & length(unique(breakpoints$breakpoint_S)) == 1) {
        breakpoints <- breakpoints %>%
          arrange(-breakpoint_S) %>%
          dplyr::first()
        breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". However S and R breakpoints are the same.")
        bp_standard <- breakpoints$site
      } else if (is.null(bp_site)) {
        breakpoints <- breakpoints %>%
          arrange(-breakpoint_S) %>%
          dplyr::first()
        breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". Using the one with the highest S breakpoint (", breakpoints$site, ").")
        bp_standard <- breakpoints$site
      } else {
        if (bp_site %in% breakpoint_sites) {
          breakpoints <- breakpoints %>% filter(site == bp_site)
          breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". Using the specified site (", bp_site, ").")
          bp_standard <- breakpoints$site
        } else {
          breakpoints <- breakpoints %>%
            arrange(-breakpoint_S) %>%
            first()
          breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". Could not find the specified site (", bp_site, "), so using the one with the highest S breakpoint (", breakpoints$site, ").")
          bp_standard <- breakpoints$site
        }
      }
    }
    breakpoint_S <- breakpoints$breakpoint_S
    breakpoint_R <- breakpoints$breakpoint_R
    if (is.na(breakpoint_S) | is.na(breakpoint_R)) {
      stop(paste("Could not determine", assay, "breakpoints using AMR package, please provide your own breakpoints"))
    }
    if (assay == "MIC") {
      breakpoint_message <- paste("MIC breakpoints determined using AMR package: S <=", breakpoint_S, "and R >", breakpoint_R)
    } else {
      breakpoint_message <- paste("Disk diffusion breakpoints determined using AMR package: S >=", breakpoint_S, "and R <", breakpoint_R)
    }
    cat(paste0("  ", breakpoint_message, "\n"))
    if (!is.na(breakpoint_message_multibp)) {
      cat(paste0("  ", breakpoint_message_multibp, "\n"))
    }
  }
  return(list(breakpoint_S = breakpoint_S, breakpoint_R = breakpoint_R, bp_standard = bp_standard))
}
