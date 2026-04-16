# subset of NCBI phenotype data that has already been re-interpreted using import_ncbi_pheno
# provided for testing solo_ppv and amr_upset functions
ecoli_pheno <- read_tsv("ecoli_pheno.tsv.gz") %>%
  import_ncbi_pheno(interpret_clsi = TRUE, interpret_ecoff = TRUE) %>%
  filter(`Scientific name` == "Escherichia coli") %>%
  filter(pheno_clsi != "NI")

# small dataframe that mimics a raw import of NCBI phenotype tab-delim file
# provided for testing import_ncbi_pheno, including re-interpreting (so only 10 rows)
ecoli_pheno_raw <- ecoli_pheno %>%
  select(-c(drug:spp_pheno)) %>%
  rename(`#BioSample` = id) %>%
  head(n = 10)

# remove unneeded column headers from phenotype object
ecoli_pheno <- ecoli_pheno %>%
  select(id:spp_pheno)

# set of AMRFinderPlus (v3.12.8, DB 2024-01-31.1) genotype results sourced from AllTheBacteria
# for the same biosamples as ecoli_pheno
# provided for testing import_amrfp, solo_ppv and amr_upset functions
ecoli_geno_raw <- read_tsv("ecoli_geno.tsv.gz")

usethis::use_data(ecoli_pheno, internal = FALSE, overwrite = TRUE)
usethis::use_data(ecoli_pheno_raw, internal = FALSE, overwrite = TRUE)
usethis::use_data(ecoli_geno_raw, internal = FALSE, overwrite = TRUE)
