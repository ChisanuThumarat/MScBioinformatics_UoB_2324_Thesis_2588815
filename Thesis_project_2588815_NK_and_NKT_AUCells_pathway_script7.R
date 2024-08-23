#T_NK_3rd_subset <- readRDS("path to base NKg")
#T_NK_3rd_subset <- readRDS("path to base NKgchemo")
#T_NK_3rd_subset <- readRDS("path to base NKgtumor")
#T_NK_3rd_subset_case_only <- readRDS("path to base NKgchemotumor_case_only")


####################################
#Input
######################################

input_seurat <- T_NK_3rd_subset

##List:
#T_NK_3rd_subset
#T_NK_3rd_subset_case_only

##################################

DefaultAssay(input_seurat) <- "RNA"

#Join all layers first
input_seurat <- JoinLayers(input_seurat)

FindVariableFeatures(input_seurat, selection.method = "vst", nfeatures = 2000, assay="RNA")

ScaleData(input_seurat)


##################################

AUC_pathway_seurat  <- input_seurat



####################################################
###Step1 Prepare gene set first 
#########################################################

msigdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
gobp_gene_sets <- msigdb_c2[msigdb_c2$gs_subcat == "CP:REACTOME", ]

go_gene_sets <- split(gobp_gene_sets$gene_symbol, gobp_gene_sets$gs_name)

#Get all HVGs from scale.data layer
all_genes <- rownames(AUC_pathway_seurat@assays$RNA$scale.data)

# Intersection with the pathways
go_gene_sets <- lapply(go_gene_sets, function(genes) {
  genes <- intersect(genes, all_genes)
  return(genes)
})

# Check the adjusted names
print(names(go_gene_sets))

length(go_gene_sets)

# Remove empty lists from the go_gene_sets
go_gene_sets <- Filter(length, go_gene_sets)

length(go_gene_sets)


####################################
#Check data layer first
####################################

expr_matrix <- GetAssayData(AUC_pathway_seurat, layer = "data")

##AUCells stage
gene_matrix <- AUCell_buildRankings(expr_matrix)

# Calculate AUC scores
suppressWarnings({
  auc_scores <- list()
  for (gs_name in names(go_gene_sets)) {
    gs <- go_gene_sets[[gs_name]]
    if (length(gs) > 0) {
      try({
        auc_scores[[gs_name]] <- AUCell_calcAUC(gs, gene_matrix)
      }, silent = TRUE) 
    } else {
      cat("Skipping empty gene set:", gs_name, "\n")
    }
  }
})

#Need to run this to satisfy R naming object
names(auc_scores) <- make.names(names(auc_scores), unique = TRUE)

# Extract AUC scores from each pathway
auc_scores_numeric <- lapply(auc_scores, function(x) {
  as.numeric(x@assays@data[[1]][1, ])
})

print(length(auc_scores_numeric[[1]]))


# Add AUCell scores to the metadata slot in Seurat (For UMAP in the future)
for (pathway in names(auc_scores_numeric)) {
  metadata_column_name <- paste0( pathway)
  AUC_pathway_seurat@meta.data[[metadata_column_name]] <- auc_scores_numeric[[pathway]]
}


# Convert AUCell scores from metadata to a matrix
auc_scores_matrix <- sapply(paste0(names(auc_scores_numeric)), function(x) {
  AUC_pathway_seurat@meta.data[[x]]
})

rownames(auc_scores_matrix) <- rownames(AUC_pathway_seurat@meta.data)

# Creating an assay with the transposed matrix
new_auc_assay <- CreateAssayObject(counts = t(auc_scores_matrix))

# Adding the new assay to the Seurat object as AUC assay
AUC_pathway_seurat[["AUC"]] <- new_auc_assay

#Find all markers but AUC
AUC_test_markers <- FindAllMarkers(AUC_pathway_seurat,assay="AUC",  min.pct = 0.25, logfc.threshold = 0.2, only.pos = TRUE)


###############################################


#####################################
#Visualization
#####################################

# To find the top markers
top_AUC_markers <- AUC_test_markers %>%
  filter(p_val_adj < 0.01) %>%  # Adjusted p-value threshold
  group_by(cluster) %>%
  top_n(n = 1, wt = abs(avg_log2FC)) 

## Top n =1 for base, 3 for chemo, 2 for tumor and 5 for chemotumor, adjusted for visualization sake. 

#List out top pathways
top_pathways <- unique(top_AUC_markers$gene)

#Adjust to make symbol work with ggplot structure
top_pathways <- gsub("-", "_", top_pathways, fixed = TRUE)

# Add active.ident to Metadata
AUC_pathway_seurat@meta.data$active_ident <- Idents(AUC_pathway_seurat)

# Extract data for melting
df_for_melt <- AUC_pathway_seurat@meta.data[, c("active_ident", top_pathways), drop = FALSE]
long_data <- melt(df_for_melt, id.vars = "active_ident", variable.name = "variable", value.name = "AUC_Score")
long_data$AUC_Score <- as.numeric(long_data$AUC_Score)

# Calculate the mean AUC score for each pathway within each NK subpopulation
average_auc_scores <- long_data %>%
  group_by(active_ident, variable) %>%
  summarise(Mean_AUC = mean(AUC_Score), .groups = 'drop')

# Define threshold for retaining pathways based on mean AUC exceeding in any subpopulation
auc_threshold_for_pathway <- 0.2
retained_pathways <- average_auc_scores %>%
  filter(Mean_AUC > auc_threshold_for_pathway) %>%
  distinct(variable) %>%
  pull(variable)

# Additional filter: Remove pathways where the mean AUC is high across ALL subclusters
high_auc_threshold <- 0.3 # Define a higher threshold for universal high expression
pathways_to_remove <- average_auc_scores %>%
  group_by(variable) %>%
  filter(all(Mean_AUC > high_auc_threshold)) %>%
  pull(variable) %>%
  unique()

# Final list of retained pathways after both filters
final_retained_pathways <- setdiff(retained_pathways, pathways_to_remove)

# Filter the data for visualization based on final retained pathways
filtered_long_data <- long_data[long_data$variable %in% final_retained_pathways,]

# Threshold for marking a cell as "expressing (active)" the pathway based on individual AUC scores
auc_threshold_for_expression <- 0.1
filtered_long_data$Expressed <- filtered_long_data$AUC_Score > auc_threshold_for_expression

# Calculate the fraction of cells expressing each gene set in each cluster
fraction_data <- filtered_long_data %>%
  group_by(active_ident, variable) %>%
  summarise(Fraction = mean(Expressed), .groups = 'drop')

# Merge the fraction data back with the average AUC data
final_plot_data <- average_auc_scores %>%
  filter(variable %in% final_retained_pathways) %>%
  inner_join(fraction_data, by = c("active_ident", "variable"))

ggplot(final_plot_data, aes(x = active_ident, y = variable, color = Mean_AUC, size = Fraction)) +
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


