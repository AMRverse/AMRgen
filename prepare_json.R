#' Generate AMR submission JSON from an AST dataset
#'
#' Converts an antimicrobial susceptibility testing (AST) dataset into
#' JSON records formatted for AMR data submission to EBI AMR portal 
#' (https://www.ebi.ac.uk/amr/amr_submission_guide/). Each row of the input
#' dataset is converted into JSON records and printed to file.
#'
#' @param ast_dataset A data frame or tibble containing antimicrobial
#' susceptibility testing results. The dataset must contain the following
#' columns: `id`, `drug_agent`, `guideline`, `method`, `measurement`,
#' `measurement_units`, `measurement_sign`, `resistance_phenotype`,
#' and `platform`.
#'
#' @param breakpoint_version Character string specifying the breakpoint
#' version used for interpretation (e.g. `"EUCAST 2024"`).
#'
#' @param submission_account Character string specifying the Webin
#' submission account identifier (e.g. `"Webin-###"`).
#'
#' @param domain Character string specifying the domain used in the
#' submission metadata (e.g. `"self.ExampleDomain"`).
#'
#' @param output_dir Character string specifying the directory where JSON
#' files should be written.
#'
#' @return Invisibly returns `NULL`. The function prints JSON-formatted
#' AMR submission records to file.
#'
#' @details
#' The function iterates over each entry in `ast_dataset` and constructs
#' a nested JSON object describing the antimicrobial susceptibility
#' testing result. Each record contains antibiotic metadata, AST
#' standards, measurement values, and resistance phenotype information.
#'
#' JSON formatting is performed using \code{jsonlite::toJSON()} with
#' `pretty = TRUE` and `auto_unbox = TRUE` .
#'
#' @examples
#' prepare_json(
#'   ast_dataset,
#'   breakpoint_version = "EUCAST 2015",
#'   submission_account = "Webin-###",
#'   domain = "self.ExampleDomain",
#'   output_dir = "/path/to/output/"
#' )
#'
#' @importFrom jsonlite write_json
#' @export
prepare_json <- function(ast_dataset,
                         breakpoint_version,
                         submission_account,
                         domain,
                         output_dir){

  # iterate through rows in tibble
  for (entry in 1:length(ast_dataset)){

    # Generate phenotype data in list format
    amr_content <- list("antibioticName" = list(value=ab_name(ast_dataset$drug_agent[entry]),
                                                iri="null"),
                        "astStandard" = list(value=ast_dataset$guideline[entry],
                                             iri="null"),
                        "breakpointVersion" = list(value=breakpoint_version,
                                                   iri="null"),
                        "laboratoryTypingMethod"=list(value=ast_dataset$method[entry],
                                                      iri="null"),
                        "measurement" = list(value=ast_dataset$measurement[entry],
                                             iri="null"),
                        "measurementUnits" = list(value=ast_dataset$measurement_units[entry],
                                                  iri="null"),
                        "measurementSign" = list(value=ast_dataset$measurement_sign[entry],
                                                 iri="null"),
                        "resistancePhenotype" = list(value=ast_dataset$resistance_phenotype[entry],
                                                     iri="null"),
                        "platform" = list(value=ast_dataset$platform[entry],
                                          iri="null"))
    # Add accession and accoutn information
    amr_record <- list("accesssion" = ast_dataset$id[entry],
                       "data" = list(list("domain" = domain,
                                          "webinSubmissionAccountId" = submission_account,
                                          "type" = "AMR",
                                          "schema" = "null",
                                          "content" = list(amr_content))))

    # Write data to json format
    write_json(amr_record, paste0(output_dir, ast_dataset$id[entry],"_", ast_dataset$drug_agent[entry],".json"), pretty = TRUE, rows=T)
  }
}
