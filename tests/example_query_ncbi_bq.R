devtools::load_all(".")

# May need to set this here (if not set elsewhere)
#if (Sys.getenv("GOOGLE_CLOUD_PROJECT") == "") {
#  Sys.setenv(GOOGLE_CLOUD_PROJECT = "my_project")
#}

message("Running query for Klebsiella + Carbapenems")
ast_raw <- query_ncbi_bq_ast(
  taxgroup = "Klebsiella pneumoniae",
  antibiotic = "meropenem"
)
message(sprintf("Found %d rows of AST data.", nrow(ast_raw)))
ast_processed <- import_ncbi_ast(ast_raw, interpret_clsi = TRUE)

amrfp <- query_ncbi_bq_geno(taxgroup = 'Klebsiella pneumoniae', geno_subclass = 'CARBAPENEM')
message(sprintf("Found %d rows of AMRFinderPlus data.", nrow(amrfp)))

geno <- import_amrfp(amrfp, sample_col = 'biosample_acc')

soloPPV <- solo_ppv_analysis(geno, ast_processed, 
                             antibiotic = 'meropenem',
                             drug_class_list = 'Carbapenems',
                             sir_col = 'pheno_clsi_mic')
upset <- amr_upset(soloPPV$amr_binary,
                   min_set_size = 5, assay = 'mic', order = 'value', 
                   print_category_counts = T)
