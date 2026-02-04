#' Co-occuring PPV Analysis
#'
#' This function performs a Positive Predictive Value (PPV) analysis for single or multiple AMR markers associated with a given antibiotic and drug class. It calculates the PPV for single and co-occurring markers and visualizes the results using UpSet-style matrix plot. 
#'
#' @param geno_table A data frame containing genotype data.
#' @param pheno_table A data frame containing phenotype data.
#' @param antibiotic String. The specific antibiotic to analyze (e.g., "Ceftazidime").
#' @param drug_class_list Vector of strings. The drug classes to filter for (e.g., c("Cephalosporins")).
#' @param geno_sample_col String. Column name for sample IDs in genotype table.
#' @param pheno_sample_col String. Column name for sample IDs in phenotype table.
#' @param sir_col String. Column name containing S/I/R data in the phenotype table.
#' @param marker_col String. Default "marker". Name of the column identifying genes.
#' @param keep_assay_values Logical. Default TRUE. Whether to keep MIC/Disk columns.
#' @param min Integer. Default 1. Minimum number of isolates required for a profile to be plotted.
#' @param axis_label_size Integer. Font size for axis text.
#' @param pd Position dodge object for plotting.
#' @param plot_cols Named vector of colors for "R" and "NWT".
#'
#' @return A list containing:
#' \item{stats}{Data frame of calculated PPV statistics per profile.}
#' \item{plot}{The combined ggplot/patchwork object.}
#' \item{data}{The filtered data frame used for the analysis.}
#' @export
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @importFrom AMR scale_fill_sir
#' @importFrom tidyr pivot_longer unnest
#' @examples
#' \dontrun{
#' geno_table <- import_amrfp(ecoli_geno_raw, "Name")
#' head(ecoli_ast)
#' cooccuringPPV_cipro <- cooccuring_ppv_analysis(
#'   geno_table = geno_table,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno_provided"
#' )
#' }

cooccuring_ppv_analysis = function (geno_table, pheno_table, antibiotic, drug_class_list, 
                                 geno_sample_col = NULL, pheno_sample_col = NULL, sir_col = NULL, 
                                 marker_col = "marker", keep_assay_values = TRUE, min = 1, 
                                 axis_label_size = 9, pd = position_dodge(width = 0.8), 
                                 plot_cols = c(R = "IndianRed", NWT = "navy")) 
{
  if (is.null(sir_col)) stop("Please specify sir_col.")
  if (!(sir_col %in% colnames(pheno_table))) stop(paste0("Column '", sir_col, "' not found."))
  
  # --- 1. Data Preparation (Same as before) ---
  amr_binary <- get_binary_matrix(geno_table, pheno_table, 
                                  antibiotic = antibiotic, drug_class_list = drug_class_list, 
                                  geno_sample_col = geno_sample_col, pheno_sample_col = pheno_sample_col, 
                                  sir_col = sir_col, keep_assay_values = keep_assay_values, 
                                  marker_col = marker_col)
  
  metadata_cols <- c("id", "pheno", "R", "NWT", "mic", "disk")
  gene_cols <- setdiff(colnames(amr_binary), metadata_cols)
  
  # Create Genotype Profiles
  amr_binary$marker <- apply(amr_binary[, gene_cols, drop=FALSE], 1, function(x) {
    active_genes <- names(x)[x == 1]
    if (length(active_genes) == 0) return(NA) 
    paste(sort(active_genes), collapse = "+") 
  })
  
  analysis_data <- amr_binary %>% 
    filter(!is.na(marker)) %>% 
    filter(!is.na(pheno))
  
  if (nrow(analysis_data) == 0) stop("No markers found matching these parameters.")
  
  # --- 2. Calculate Stats ---
  stats_R <- analysis_data %>% group_by(marker) %>% 
    summarise(x = sum(R, na.rm = T), n = n(), p = x/n, 
              se = sqrt(p * (1 - p)/n), 
              ci.lower = max(0, p - 1.96 * se), 
              ci.upper = min(1, p + 1.96 * se)) %>% 
    mutate(category = "R")
  
  stats_NWT <- analysis_data %>% group_by(marker) %>% 
    summarise(x = sum(NWT), n = n(), p = x/n, 
              se = sqrt(p * (1 - p)/n), 
              ci.lower = max(0, p - 1.96 * se), 
              ci.upper = min(1, p + 1.96 * se)) %>% 
    mutate(category = "NWT")
  
  profile_stats <- bind_rows(stats_R, stats_NWT) %>% 
    relocate(category, .before = x) %>% rename(ppv = p)
  
  # --- 3. Prepare Plot Data ---
  
  # Filter by minimum count
  markers_to_keep <- profile_stats %>% 
    filter(category == "R") %>% # Use R stats to filter/sort
    filter(n >= min) %>%
    arrange(n) %>% # Sort by frequency (Low to High) so common ones are at the top of plot
    pull(marker)
  
  if(length(markers_to_keep) == 0) stop("No combinations pass the 'min' threshold.")
  
  # Filter main datasets
  analysis_data <- analysis_data %>% filter(marker %in% markers_to_keep)
  profile_stats <- profile_stats %>% filter(marker %in% markers_to_keep)
  
  # Ensure Factor Levels match across all plots (Crucial for alignment)
  analysis_data$marker <- factor(analysis_data$marker, levels = markers_to_keep)
  profile_stats$marker <- factor(profile_stats$marker, levels = markers_to_keep)
  
  # --- 4. Build The Matrix Plot (The "UpSet" dots) ---
  
  # Deconstruct the "GeneA+GeneB" strings back into data for dots
  matrix_data <- data.frame(marker = markers_to_keep) %>% 
    mutate(marker = factor(marker, levels = markers_to_keep)) %>%
    mutate(genes_list = strsplit(as.character(marker), "\\+")) %>%
    unnest(genes_list) %>%
    rename(gene = genes_list)
  
  # Sort Genes on X-axis by prevalence
  gene_order <- matrix_data %>% count(gene) %>% arrange(desc(n)) %>% pull(gene)
  matrix_data$gene <- factor(matrix_data$gene, levels = gene_order)
  
  # Calculate line segments (connect min gene to max gene for each row)
  matrix_lines <- matrix_data %>% 
    group_by(marker) %>% 
    summarise(min = min(as.numeric(gene)), max = max(as.numeric(gene)))
  
  g_matrix <- ggplot(matrix_data, aes(x = gene, y = marker)) + 
    # The Matrix Grid
    geom_segment(data = matrix_lines, aes(x = min, xend = max, y = marker, yend = marker), 
                 color = "grey50", linewidth = 0.5) +
    geom_point(size = 3) + 
    theme_light() + 
    scale_x_discrete(position = "top") + # Put gene names at top
    labs(x = "", y = "") +
    theme(axis.text.y = element_blank(), # Hide marker text (the dots ARE the text now)
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 0), # Rotate gene names
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  # --- 5. Build Phenotype Plot (Middle) ---
  g_pheno <- ggplot(analysis_data, aes(y = marker, fill = pheno)) + 
    geom_bar(stat = "count", position = "fill") + 
    scale_fill_sir() + 
    geom_text(aes(label = after_stat(count)), stat = "count", 
              position = position_fill(vjust = 0.5), size = 3) + 
    theme_light() + 
    labs(y = "", x = "Proportion", fill = "Phenotype") +
    theme(axis.text.y = element_blank(), # Hide Y labels (Matrix handles this)
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = axis_label_size),
          legend.position = "bottom")
  
  # --- 6. Build PPV Plot (Right) ---
  
  # Corrected label logic 
  # 1. We need a lookup table that maps each marker to its 'n' (count)
  # 2. We explicitly arrange it by the factor levels used in the plot
  label_data <- profile_stats %>%
    filter(category == "R") %>%
    mutate(marker = factor(marker, levels = markers_to_keep)) %>% # Enforce same order
    arrange(marker) # Sort by that order
  
  # Create the label vector in the correct order
  y_axis_labels <- paste0("(n=", label_data$n, ")")
  
  g_ppv <- ggplot(profile_stats, aes(y = marker, group = category, col = category)) + 
    geom_vline(xintercept = 0.5, linetype = 2, color = "grey") + 
    geom_linerange(aes(xmin = ci.lower, xmax = ci.upper), position = pd) + 
    geom_point(aes(x = ppv), position = pd) + 
    theme_light() + 
    # Use the ordered becypr here:
    scale_y_discrete(labels = y_axis_labels, position = "right") + 
    labs(y = "", x = "PPV", col = "Category") + 
    scale_colour_manual(values = plot_cols) + 
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    theme(axis.text.x = element_text(size = axis_label_size),
          axis.text.y = element_text(size = axis_label_size), 
          legend.position = "bottom")
  
  # --- 7. Assembly ---
  header <- paste("Genotype Profiles:", paste0(drug_class_list, collapse = ", "))
  
  # Combine with patchwork
  # Layout: Matrix | Pheno | PPV
  combined_plot <- (g_matrix + g_pheno + g_ppv) + 
    plot_layout(widths = c(1.2, 1, 1), guides = "collect") + 
    plot_annotation(title = header, subtitle = paste("vs phenotype for:", antibiotic)) & 
    theme(legend.position = "bottom")
  
  print(combined_plot)
  return(list(stats = profile_stats, plot = combined_plot, data = analysis_data))
}