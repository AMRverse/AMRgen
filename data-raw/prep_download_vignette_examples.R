# Main vignette
ecoli_cip_mic_data <- get_eucast_mic_distribution("cipro", "E. coli")
usethis::use_data(ecoli_cip_mic_data, internal = FALSE, overwrite = TRUE)

ecoli_cip <- ecoli_ast$mic[ecoli_ast$drug_agent == "CIP"]
ecoli_cip_vs_ref <- compare_mic_with_eucast(ecoli_cip, ab = "cipro", mo = "E. coli")
usethis::use_data(ecoli_cip_vs_ref, internal = FALSE, overwrite = TRUE)


# Download data vignette
staph_ast_ncbi <- download_ncbi_ast(
  species = "Staphylococcus aureus",
  antibiotic = c("amikacin", "DOX"), # antibiotics can be listed in short or long form
  reformat = TRUE,
  interpret_eucast = TRUE
) # reformat must be true to use interpret_* argument
usethis::use_data(staph_ast_ncbi, internal = FALSE, overwrite = TRUE)

staph_ast_ncbi_raw <- download_ncbi_ast(
  species = "Staphylococcus aureus",
  antibiotic = c("amikacin", "DOX"),
  reformat = FALSE,
  interpret_eucast = FALSE
)
usethis::use_data(staph_ast_ncbi_raw, internal = FALSE, overwrite = TRUE)


# downloading via google bigquery - note you may need to add project_id to run
staph_ast_ncbi_cloud_raw <- query_ncbi_bq_ast(
  taxgroup = "Staphylococcus aureus",
  antibiotic = c("amikacin", "DOX")
)
usethis::use_data(staph_ast_ncbi_cloud_raw, internal = FALSE, overwrite = TRUE)

staph_geno_ncbi_cloud_raw <- query_ncbi_bq_geno(
  taxgroup = "Staphylococcus aureus",
  geno_class = c("AMINOGLYCOSIDE", "TETRACYCLINE")
)
# file is too big including all samples that have ANY AST data in NCBI, narrow down to those with
staph_geno_ncbi_cloud_raw <- staph_geno_ncbi_cloud_raw %>% filter(biosample_acc %in% staph_ast_ncbi_cloud_raw$BioSample)
usethis::use_data(staph_geno_ncbi_cloud_raw, internal = FALSE, overwrite = TRUE)


staph_ast_ebi <- download_ebi(
  genus = "Staphylococcus",
  antibiotic = c("amikacin", "DOX"),
  reformat = TRUE,
  interpret_eucast = TRUE,
  interpret_clsi = TRUE,
  interpret_ecoff = TRUE
) # chose which guideline to use for re-interpretation. reformat must be TRUE for re-interpretation
usethis::use_data(staph_ast_ebi, internal = FALSE, overwrite = TRUE)

staph_geno_ebi <- download_ebi(
  data = "genotype",
  genus = "Staphylococcus",
  geno_class = c("AMINOGLYCOSIDE", "TETRACYCLINE"),
  reformat = T
)
usethis::use_data(staph_geno_ebi, internal = FALSE, overwrite = TRUE)
