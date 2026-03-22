kp_mero <- download_ebi(
  antibiotic = "meropenem",
  species = "Klebsiella pneumoniae",
  reformat = TRUE,
  interpret_eucast = TRUE,
  interpret_ecoff = TRUE
)

# Filter for isolates in EuSCAPE paper (PMID: 31358985)
kp_mero_euscape <- kp_mero %>% filter(grepl("31358985", source))

# There are assemblies from NCBI that are flagged for contamination and supposed to be excluded. For example, see SAMEA3729690 (https://www.ncbi.nlm.nih.gov/datasets/genome/?biosample=SAMEA3729690)

contaminated_assemblies <- c("SAMEA3729690", "SAMEA3721062", "SAMEA3721052", "SAMEA3720966", "SAMEA3673128", "SAMEA3538742", "SAMEA3721188", "SAMEA3649589", "SAMEA3538652", "SAMEA3649503", "SAMEA3538911", "SAMEA3727711", "SAMEA3649452", "SAMEA3649453", "SAMEA3649454", "SAMEA3649467", "SAMEA3721063", "SAMEA3538862", "SAMEA3538667", "SAMEA3673004", "SAMEA3729818", "SAMEA3729660", "SAMEA3673078", "SAMEA3673097")

# Remove contaminated assemblies from phenotype list
kp_mero_euscape <- kp_mero_euscape %>%
  filter(!id %in% contaminated_assemblies)

usethis::use_data(kp_mero_euscape, internal = FALSE, overwrite = TRUE)


kp_mero_amrfp <- download_ebi(
  data = "genotype", species = "Klebsiella pneumoniae",
  reformat = T
)

# Filter for isolates in EuSCAPE paper with meropenem phenotypes and remove contaminated samples
kp_mero_amrfp <- kp_mero_amrfp %>% filter(id %in% kp_mero_euscape$id)

# There are assemblies from NCBI that are supposed to be excluded and flagged for contamination. For example, see SAMEA3729690 (https://www.ncbi.nlm.nih.gov/datasets/genome/?biosample=SAMEA3729690)

kp_mero_amrfp <- kp_mero_amrfp %>%
  filter(!id %in% contaminated_assemblies)

usethis::use_data(kp_mero_amrfp, internal = FALSE, overwrite = TRUE)