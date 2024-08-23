#########################################
#AUCells-based on fucntional markers
#########################################

#T_NK_3rd_subset <- readRDS("path to base NKg")
#T_NK_3rd_subset <- readRDS("path to base NKgchemo")
#T_NK_3rd_subset <- readRDS("path to base NKgtumor")
#T_NK_3rd_subset_case_only <- readRDS("path to base NKgchemotumor_case_only")

# Define the gene sets
cytotoxicity_markers <- c("GNLY", "GZMA", "GIMAP7", "PSMB10", "TBC1D10C", "FGFBP2", "ACTB", "CSK")
activation_markers <- c("CD69", "NCR1", "KLRK1", "CD226", "CD96", "FOS", "IFNG", "ATF3", "PELI2", "IGF1R", "IL12RB2", "IER3", "EGR1", "AREG")
adhesion_and_migration_markers <- c("CD44", "SELL", "CXCR6", "CHPT1", "TIAM1")
proliferation_markers <- c("MKI67", "PVT1", "PBX4", "IRS2", "IL7R")
inhibitory_markers <- c("KLRD1", "KLRC1", "KLRC2", "PDCD1", "HAVCR2", "KIR2DL4", "BACH2", "NFKBID", "IKZF2")
chemokines_and_receptors_markers <- c("XCL1", "XCL2", "CCL3", "CCL4", "CCL4L2", "CCL3L1")
stress_and_exhaustion_markers <- c("HSPA1A", "HSP90AA1", "DNAJB1", "HSPA1B", "HSPH1", "HSPD1", "FKBP4", "HSPA6", "HSPB1", "ZFAND2A")
general_nk_cell_markers <- c("CD7", "NCAM1", "FCGR3A", "B3GAT1", "KLRF1", "NFKB1", "P2RY8", "PHLDB2", "BNC2", "AXIN1", "INPP5A")
miscellaneous_other_functional_roles <- c("CRTAM", "LINC02446", "RILPL2", "MAML3", "B3GNT7", "SSBP2", "ABTB2", "FTH1", "C1orf21", "TGFBR3", "KLF3", "PEX14", "OASL", "WDR47", "VPS13D")


# Combine all gene sets into a list
gene_sets <- list(
  "Cytotoxicity" = cytotoxicity_markers,
  "Activation" = activation_markers,
  "Adhesion_and_Migration" = adhesion_and_migration_markers,
  "Proliferation" = proliferation_markers,
  "Inhibitory" = inhibitory_markers,
  "Chemokines_and_Receptors" = chemokines_and_receptors_markers,
  "Stress_and_Exhaustion" = stress_and_exhaustion_markers,
  "General_NK_Cell" = general_nk_cell_markers,
  "Miscellaneous_Other_Functional_Roles" = miscellaneous_other_functional_roles
)

#Obtain expression matrix
expr_matrix <- as.matrix(GetAssayData(T_NK_3rd_subset, layer = "data"))



##AUCells steps
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)
auc_matrix <- as.matrix(getAUC(cells_AUC))

# Add AUC scores to the metadata of Seurat
T_NK_3rd_subset <- AddMetaData(T_NK_3rd_subset, metadata = t(auc_matrix))

# Plot UMAP colored by AUC scores of each gene set
for (gene_set in names(gene_sets)) {
  umap_plot <- FeaturePlot(T_NK_3rd_subset, features = gene_set) + 
    scale_color_viridis_c() + 
    ggtitle(paste("UMAP of", gene_set, "AUC Scores"))
  
  print(umap_plot)
}

####################################

top_pathways <- unique(names(gene_sets))

# Add active.ident to Metadata
T_NK_3rd_subset@meta.data$active_ident <- Idents(T_NK_3rd_subset)

# Extract data for melting
df_for_melt <- T_NK_3rd_subset@meta.data[, c("active_ident", top_pathways), drop = FALSE]
long_data <- melt(df_for_melt, id.vars = "active_ident", variable.name = "variable", value.name = "AUC_Score")
long_data$AUC_Score <- as.numeric(long_data$AUC_Score)

# Calculate the mean AUC score 
average_auc_scores <- long_data %>%
  group_by(active_ident, variable) %>%
  summarise(Mean_AUC = mean(AUC_Score), .groups = 'drop')

# Threshold for retaining pathways based on mean AUC exceeding in any subpopulation
auc_threshold_for_pathway <- 0
retained_pathways <- average_auc_scores %>%
  filter(Mean_AUC > auc_threshold_for_pathway) %>%
  distinct(variable) %>%
  pull(variable)

# Additional filter: Remove pathways where the mean AUC is high across ALL subclusters
high_auc_threshold <- 1 # Define a higher threshold for universal high expression
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
auc_threshold_for_expression <- 0
filtered_long_data$Expressed <- filtered_long_data$AUC_Score > auc_threshold_for_expression

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

#####################################
