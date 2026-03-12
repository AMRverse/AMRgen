# subset of AST data that has already been re-interpreted using import_ncbi_ast
ast_CLI_public <- read_tsv("data-raw/ast_CLI_public.tsv.gz")

# provide AMRFinderPlus results
afp_CLI_public <- read_tsv("data-raw/afp_CLI_public.tsv.gz")

# provide ST data from AllTheBacteria
ST_data_CLI <- read_tsv("data-raw/ST_data_CLI.tsv.gz")


usethis::use_data(ast_CLI_public, internal = FALSE, overwrite = TRUE)
usethis::use_data(afp_CLI_public, internal = FALSE, overwrite = TRUE)
usethis::use_data(ST_data_CLI, internal = FALSE, overwrite = TRUE)
