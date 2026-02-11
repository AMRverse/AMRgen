amrfp_drugs_table <- read_tsv("amrfp_drug_classes_agents.tsv")

usethis::use_data(amrfp_drugs_table, internal = FALSE, overwrite = TRUE)
