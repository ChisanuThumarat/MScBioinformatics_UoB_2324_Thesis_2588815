#T_NK_3rd_subset <- readRDS("path to base NKg")
#T_NK_3rd_subset <- readRDS("path to base NKgchemo")
#T_NK_3rd_subset <- readRDS("path to base NKgtumor")
#T_NK_3rd_subset_case_only <- readRDS("path to base NKgchemotumor_case_only")

AUC_regulon_seurat <- T_NK_3rd_subset 
#T_NK_3rd_subset_case_only

AUC_regulon_seurat <- JoinLayers(AUC_regulon_seurat)
#No need to rerun scale.data due to not using HVGs


######################################
#Regulons AUCells
#####################################

#Preparing C3 pathway list
msigdb_c3 <- msigdbr(species = "Homo sapiens", category = "C3")
msigdb_c3_TF_Legacy <- msigdb_c3[msigdb_c3$gs_subcat == "TFT:TFT_Legacy", ]
msigdb_c3_GTRD <- msigdb_c3[msigdb_c3$gs_subcat == "TFT:GTRD", ]
msigdb_c3_TF <- rbind(msigdb_c3_TF_Legacy, msigdb_c3_GTRD)

# Prepare the gene sets
c3_gene_sets <- msigdb_c3_TF$gene_symbol
regulon_list <- unique(c(c3_gene_sets))

####### Preparation info step #########
all_genes_data <- rownames(AUC_regulon_seurat@assays$RNA$data)
regulon_list <- intersect(regulon_list, all_genes_data)



################################################
##A bit reuse the word due to originally create it via pathway version so "go_regulon_sets = regulon_sets" alone

go_regulon_sets <- split(msigdb_c3_TF$gene_symbol, msigdb_c3_TF$gs_name)

#Intersection
go_regulon_sets <- lapply(go_regulon_sets, function(genes) {
  genes <- intersect(genes, regulon_list)  # Keep only genes that are part of the regulons
  return(genes)
})

# Get all genes from gene matrix data
all_genes_data <- rownames(AUC_regulon_seurat@assays$RNA$data)

go_regulon_sets <- lapply(go_regulon_sets, function(genes) {
  genes <- intersect(genes, all_genes_data)
  return(genes)
})

# Remove empty lists from the go_regulon_sets
go_regulon_sets <- Filter(length, go_regulon_sets)

# Extracting the regulon-supported genes from C3
regulon_supported_genes <- unique(unlist(go_regulon_sets))

# Filter the expression matrix to only include these genes
expr_matrix2 <- GetAssayData(AUC_regulon_seurat, layer = "data")
filtered_expr_matrix <- expr_matrix2[regulon_supported_genes, , drop = FALSE]

###################AUC stage

gene_matrix2 <- AUCell_buildRankings(filtered_expr_matrix)
suppressWarnings({
  auc_scores2 <- list()
  for (gs_name in names(go_regulon_sets)) {
    gs <- go_regulon_sets[[gs_name]]
    if (length(gs) > 0) {
      auc_scores2[[gs_name]] <- AUCell_calcAUC(gs, gene_matrix2)
    } else {
      cat("Skipping empty gene set:", gs_name, "\n")
    }
  }
})

names(auc_scores2) <- make.names(names(auc_scores2), unique = TRUE)

auc_scores_numeric2 <- lapply(auc_scores2, function(x) {
  as.numeric(x@assays@data[[1]][1, ])
})


print(length(auc_scores_numeric2[[1]]))


# Add AUCell scores directly to the metadata in Seurat
for (pathway in names(auc_scores_numeric2)) {
  metadata_column_name <- paste0( pathway)
  AUC_regulon_seurat@meta.data[[metadata_column_name]] <- auc_scores_numeric2[[pathway]]
}

# Convert AUCell scores from metadata to a matrix
auc_scores_matrix2 <- sapply(paste0(names(auc_scores_numeric2)), function(x) {
  AUC_regulon_seurat@meta.data[[x]]
})

rownames(auc_scores_matrix2) <- rownames(AUC_regulon_seurat@meta.data)

# Creating an assay with the transposed matrix
new_auc_assay2 <- CreateAssayObject(counts = t(auc_scores_matrix2))

# Adding the AUC assay into Seurat
AUC_regulon_seurat[["AUC"]] <- new_auc_assay2

save(list = ls(), file = "AUCells_NKg_base_regulon.RData")

#Find all markers but AUC
AUC_regulon_markers <- FindAllMarkers(AUC_regulon_seurat,assay="AUC", min.pct = 0.25, logfc.threshold = 0.2, only.pos = TRUE)

# Save all objects in the environment
save(list = ls(), file = "AUCells_NKg_base_regulon.RData")



#####################



#####################################
#Visualization - regulon
#####################################

# To find the top regulons
top_AUC_markers2 <- AUC_regulon_markers %>%
  filter(p_val_adj < 0.01) %>%  
  group_by(cluster) %>%
  top_n(n = 2, wt = abs(avg_log2FC)) 

## top number adjusted for visualization sake. 

top_pathways2 <- unique(top_AUC_markers2$gene)
top_pathways2 <- gsub("-", "_", top_pathways2, fixed = TRUE)
AUC_regulon_seurat@meta.data$active_ident <- Idents(AUC_regulon_seurat)

# Extract data for melting
df_for_melt <- AUC_regulon_seurat@meta.data[, c("active_ident", top_pathways2), drop = FALSE]
long_data <- melt(df_for_melt, id.vars = "active_ident", variable.name = "variable", value.name = "AUC_Score")
long_data$AUC_Score <- as.numeric(long_data$AUC_Score)

# Calculate the mean AUC score for each targeted motifs within each NK subpopulation
average_auc_scores <- long_data %>%
  group_by(active_ident, variable) %>%
  summarise(Mean_AUC = mean(AUC_Score), .groups = 'drop')

# Define threshold for retaining targeted motifs based on mean AUC exceeding in any subpopulation
auc_threshold_for_pathway <- 0.1
retained_pathways <- average_auc_scores %>%
  filter(Mean_AUC > auc_threshold_for_pathway) %>%
  distinct(variable) %>%
  pull(variable)

# Additional filter: Remove pathways where the mean AUC is high across ALL subclusters
high_auc_threshold <- 0.2 # Define a higher threshold for universal high expression
pathways_to_remove <- average_auc_scores %>%
  group_by(variable) %>%
  filter(all(Mean_AUC > high_auc_threshold)) %>%
  pull(variable) %>%
  unique()

# Final list of retained pathways after both filters
final_retained_pathways <- setdiff(retained_pathways, pathways_to_remove)

# Filter the data for visualization based on final retained pathways
filtered_long_data <- long_data[long_data$variable %in% final_retained_pathways,]

# Threshold for marking a cell as "expressing" the pathway based on individual AUC scores
auc_threshold_for_expression <- 0.1
filtered_long_data$Expressed <- filtered_long_data$AUC_Score > auc_threshold_for_expression

fraction_data <- filtered_long_data %>%
  group_by(active_ident, variable) %>%
  summarise(Fraction = mean(Expressed), .groups = 'drop')

# Merge the fraction data back with the average AUC data
final_plot_data2 <- average_auc_scores %>%
  filter(variable %in% final_retained_pathways) %>%
  inner_join(fraction_data, by = c("active_ident", "variable"))

ggplot(final_plot_data2, aes(x = active_ident, y = variable, color = Mean_AUC, size = Fraction)) +
  geom_point() +
  scale_color_viridis_c() +  
  scale_size(range = c(2, 12)) +
  labs(title = "Top Pathway Enrichment Scores by NK Cell Cluster",
       x = "NK Cell Cluster",
       y = "Pathway",
       color = "Mean AUC Score",
       size = "Fraction of Cells Expressing") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
###############################################








