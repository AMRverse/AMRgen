# subset of phenotype data that has already been re-interpreted using import_ncbi_pheno
pheno_CLI_public <- format_pheno("data-raw/ast_CLI_public.tsv.gz", ab_col = "drug_agent")

# provide AMRFinderPlus results
afp_CLI_public <- read_tsv("data-raw/afp_CLI_public.tsv.gz") %>%
  rename(id = Name)

# provide ST data from AllTheBacteria
ST_data_CLI <- read_tsv("data-raw/ST_data_CLI.tsv.gz") %>%
  rename(id = Name)


usethis::use_data(pheno_CLI_public, internal = FALSE, overwrite = TRUE)
usethis::use_data(afp_CLI_public, internal = FALSE, overwrite = TRUE)
usethis::use_data(ST_data_CLI, internal = FALSE, overwrite = TRUE)
