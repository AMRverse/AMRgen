# E. coli data for the Concordance vignette from Mills et al, Genome Medicine (2022) 14:147

library(readr)
# Phenotypic data obtained from the EBI AMR portal with download_ebi function, reformatted but NOT reinterpreted
pheno_eco_2075 <- read_tsv("/Users/silargi/Documents/ESGEM-AMR/AMRgen/data-raw/pheno_eco_2075.tsv.gz") 

# Phenotypic data obtained from All The Bacteria
geno_eco_2075 <- read_tsv("/Users/silargi/Documents/ESGEM-AMR/AMRgen/data-raw/geno_eco_2075.tsv.gz")

usethis::use_data(pheno_eco_2075, internal = FALSE, overwrite = TRUE)
usethis::use_data(geno_eco_2075, internal = FALSE, overwrite = TRUE)
