# Load raw geno/pheno file
salm_raw <- read_csv("data-raw/Salmonella_pheno_geno_data.csv")

usethis::use_data(salm_raw, internal = FALSE, overwrite = TRUE)
