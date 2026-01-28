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

# # Download EBI data
# # Zoe A. Dyson (zoe.dyson@lshtm.ac.uk)
# # Last updated 28/01/2026
#
# # List required packages
# packages <- c("arrow", "RCurl)
#
# # Install missing packages
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }
#
# # Load packages
# invisible(lapply(packages, library, character.only = TRUE))

#' @title import_ebi
#'
#' @param user_genus
#' @param user_release
#' @param user_antibiotic_name
#'
#' @return A data frame containing EBI genotype data
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' import_ebi(
#'     user_genus="Salmonella",
#'     user_release="2025-12",
#'     user_antibiotic_name="ampicillin"
#' )
#' }
import_ebi <- function(user_genus=NULL, user_release=NULL, user_antibiotic_name=NULL) {
  
  # EBI source url
  ebi_url <- "ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/"

  if (!is.null(user_release)){
    ebi_geno <- getBinaryURL(print(gsub("ftp\\:", "https\\:", paste0(ebi_url,user_release,"/","genotype.parquet"))))
    ebi_geno <- read_parquet(ebi_geno)
    
  }else{
  
    # get list of releases
    folders <- str_split(getURL(ebi_url, dirlistonly = TRUE), "\n")[[1]]
    # get latest
    latest_release <- folders[!folders %in% c("releases.yml", "")] %>% 
      max()
    ebi_geno <- getBinaryURL(gsub("ftp\\:", "https\\:", paste0(ebi_url,latest_release,"/","genotype.parquet")))
    ebi_geno <- read_parquet(ebi_geno)
    
  }
  
  # filter by genus if required
  if (!is.null(user_genus)){
    ebi_geno <- ebi_geno %>% 
      filter(genus==user_genus)
  }
  
  # filter by genus if required
  if (!is.null(user_antibiotic_name)){
    ebi_geno <- ebi_geno %>% 
      filter(antibiotic_name==user_antibiotic_name)
  }
  
  # return data frame
  return(ebi_geno)
}
