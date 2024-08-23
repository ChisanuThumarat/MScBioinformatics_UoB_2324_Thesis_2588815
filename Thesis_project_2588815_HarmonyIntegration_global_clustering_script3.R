
#Set seed
set.seed(42) 

###########################################
#Post-merging and doubletfinder
###########################################

##Update seurat phenotype before analysis

# Create the phenotype data frame
phenotype_data <- data.frame(
  sample_ID = c('s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10',
                's11', 's12', 's13', 's14', 's15', 's16', 's17', 's18', 's19', 's20',
                's21', 's22', 's23', 's24', 's25'),
  chemo = c('CT', 'TN', 'CT', 'TN', 'CT', 'CT', 'TN', 'CT', 'TN', 'TN',
            'CT', 'TN', 'CT', 'TN', 'TN', 'CT', 'TN', 'TN', 'CT', 'CT',
            'TN', 'TN', 'TN', 'TN', 'TN'),
  tumorType = c('HGSOC', 'HGSOC', 'HGSOC', 'GI_MET', 'Norm', 'HGSOC', 'HGSOC', 'HGSOC', 'HGSOC', 'HGSOC',
                'EAC', 'EAC', 'EAC', 'EAC', 'EAC', 'EAC', 'Norm', 'EAC', 'Norm', 'EAC',
                'PDAC', 'Benign', 'PDAC', 'PDAC', 'PDAC')
)


metadata <- combined_seurat@meta.data

metadata <- metadata %>%
  left_join(phenotype_data, by = c("sample_ID" = "sample_ID"))

# Update the Seurat  
combined_seurat <- AddMetaData(combined_seurat, metadata)


#########################################
#Harmony integration and SCTransform
########################################

#FindVariableFeatures first for top variable features and PCA
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000, assay="RNA")
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000, assay="SCT")

#Mitochondrial percentage for SCT regression 
combined_seurat <- PercentageFeatureSet(combined_seurat, pattern = "^MT-", col.name = "percent.mt")

#SCTransform normalization
combined_seurat <- SCTransform(combined_seurat, vars.to.regress = "percent.mt", verbose = FALSE)

DefaultAssay(combined_seurat) <- "SCT"

combined_seurat <- RunPCA(combined_seurat, assay.use = "SCT", verbose = FALSE)

#Rerun Harmony for batch correction
combined_seurat <- RunHarmony(combined_seurat, assay.use = "SCT",  "sample_ID")

#Join all 25 Seurat layers
combined_seurat <- JoinLayers(combined_seurat)

# Save the integrated Seurat object
#saveRDS(combined_seurat, "SCT_integrated_seurat.rds")

#######################
# Load the data
#combined_seurat <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/SCT_integrated_seurat.rds")

#Meta_data <- read.csv("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/DATA/metadata_cr_v7.csv")
#######################


########################################################
#FindNeighbors and FindClusters + UMAP
########################################################

DefaultAssay(combined_seurat) <- "SCT"

# Run UMAP 
combined_seurat <- RunUMAP(combined_seurat, assay = "SCT",  reduction = "harmony", dims = 1:50)

# Find Neighbors and Clusters
combined_seurat <- FindNeighbors(combined_seurat, assay = "SCT", reduction = "harmony", dims = 1:50)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

DimPlot(combined_seurat, reduction = "umap", group.by = "sample_ID") + ggtitle("UMAP Plot Combined Seurat")
DimPlot(combined_seurat, reduction = "harmony", group.by = "sample_ID") + ggtitle("Harmony Plot Combined Seurat")

# Save the integrated Seurat object
#saveRDS(combined_seurat, "join_layer_integrated_seurat.rds")


######################################
# Load the data
#combined_seurat <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/join_layer_integrated_seurat.rds")

#Meta_data <- read.csv("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/DATA/metadata_cr_v7.csv")
########################################

#######################################################
#Unsupervised clustering and clustering defining
#######################################################

###Major Lineage checking

#Dotplot
MajorTypes_panel=c("CD3D",# T cells
                   "MS4A1", # B cells
                   "IGKC", # Plasmablasts
                   "EPCAM", # Epithelial
                   "MKI67", # Cycling
                   "PECAM1", # Endothelial
                   "DCN", # Fibroblasts
                   "COL1A1", # Fibroblast and endothelial
                   "PDPN", # Fibroblasts
                   "LYZ", # Myeloid
                   "TPSAB1", # Mast
                   "GNLY", # NK and T
                   "XCL1", # NK
                   "FGFBP2") # NK
DotPlot(object = combined_seurat, features = MajorTypes_panel) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")




#Outputs

#"CD3D"# T cells  : 0, 1, 3, 5, 6, 7, 8, 13, 18, 21, 24, 25, 26, 27
#"MS4A1", # B cells : 4, 17, 19, 28, 29
#"IGKC", # Plasmablasts : 4, 17, 19
#"EPCAM", # Epithelial : 16
#"MKI67", # Cycling : 22
#"PECAM1", # Endothelial : 9, 10, 11, 12, 15, 16, 22, 23
#"DCN", # Fibroblasts :  2, 10, 16, 18, 25
#"COL1A1", # Fibroblast and endothelial : 2, 10, 16, 18, 25
#"PDPN", # Fibroblasts : 2, 10, 16, 18, 25
#"LYZ", # Myeloid : 9, 11, 14, 15, 23
#"TPSAB1", # Mast : 20
#"GNLY", # NK and T : 0, 8, 13
#"XCL1", # NK :0, 8, 13
#"FGFBP2"# NK :0, 8, 13


###LYZ and others checking
checking_panels=c("LYZ", "CD14", "FCGR3A", "NKG7") 
DotPlot(object = combined_seurat, features = checking_panels)  + scale_colour_gradient2(low = "blue", mid = "white", high = "red")


###General NK panel checking
nk_markers <- c("NCAM1", "KLRK1" , "KLRF1", "NKG7", "GNLY", "FCGR3A", "CD3D", "CD3G", "CD3E")
DotPlot(object = combined_seurat, features = nk_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")


#"MS4A1", # B cells : 4, 12, 14, 17
#"IGKC", # Plasmablasts : 4, 17, 19, 28, 29
#"EPCAM", # Epithelial : 16
#"MKI67", # Cycling : 22
#"PECAM1", # Endothelial : 9, 10, 11, 12, 15, 16, 22, 23
#"DCN", # Fibroblasts :  2, 10, 16, 18, 25
#"COL1A1", # Fibroblast and endothelial : 2, 10, 16, 18, 25
#"PDPN", # Fibroblasts : 2, 10, 16, 18, 25
#"LYZ", # Myeloid : 9, 11, 14, 15, 23
#"TPSAB1", # Mast : 20

#T_cells  :  1, 3, 5, 6, 7, 18, 21, 24, 25, 26, 27
#B cells : 4
#Plasmablasts : 17, 19, 28, 29
#Epithelial : 16
#Cycling : 22
#Endothelial :  10
#Fibroblasts :  2, 18, 25
#Myeloid : 9, 11, 14, 15, 23, 12
#Mast : 20
#NK: 8, 13
#NKT: 0

###Defining clusters
new.cluster.ids <- c(
  "NKT",  # 0
  "T_cells",  # 1
  "Fibroblasts",  # 2
  "T_cells",  # 3
  "B_cells",  # 4
  "T_cells",  # 5
  "T_cells",  # 6
  "T_cells",  # 7
  "NK",  # 8
  "Myeloid",  # 9
  "Endothelial",  # 10
  "Myeloid",  # 11
  "Myeloid",  # 12
  "NK",  # 13
  "Myeloid",  # 14
  "Myeloid",  # 15
  "Epithelial",  # 16
  "Plasmablasts",  # 17
  "T_cells",  # 18
  "Plasmablasts",  # 19
  "Mast",  # 20
  "T_cells",  # 21
  "Cycling",  # 22
  "Myeloid",  # 23
  "T_cells",  # 24
  "Fibroblasts",  # 25
  "T_cells",  # 26
  "T_cells",  # 27
  "Plasmablasts",  # 28
  "Plasmablasts"  # 29
)

# Rename the clusters 
names(new.cluster.ids) <- levels(combined_seurat)
combined_seurat <- RenameIdents(combined_seurat, new.cluster.ids)

# Plot the UMAP check
DimPlot(combined_seurat, reduction = "umap", label = TRUE, pt.size = 0.5)



# Save the integrated Seurat object
#saveRDS(combined_seurat, "Assigned_label_first_round_integrated_seurat.rds")


###Subset
T_NK_cell_subset <- subset(combined_seurat, idents = c("NK", "NKT"))

# Save 
saveRDS(T_NK_cell_subset, "T_NK_cell_subset.rds")

##############################################
# Load the data
#combined_seurat <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/Assigned_label_first_round_integrated_seurat.rds")
#Set back to SCT
#DefaultAssay(combined_seurat) <- "SCT"
##############################################


#############################
#Visualization for CT and TN
############################

color_palette_ct <- c("CT" = "salmon", "TN" = "white")

DimPlot(combined_seurat, reduction = "umap", group.by = "chemo.y", cols = color_palette_ct) + 
  ggtitle("UMAP highlighting CT cells")

color_palette_tn <- c("CT" = "white", "TN" = "cyan")

DimPlot(combined_seurat, reduction = "umap", group.by = "chemo.y", cols = color_palette_tn) + 
  ggtitle("UMAP highlighting CT cells")

#############################################
#Major lineage and NK markers visualization
#############################################

general_nk_cell_markers <- c("CD7", "NCAM1", "FCGR3A","B3GAT1", "KLRF1", "NFKB1", "P2RY8",  "AXIN1")
DotPlot(object = combined_seurat, features = general_nk_cell_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+ RotatedAxis()+coord_flip()

nk_markers <- c("NCAM1", "KLRK1" , "KLRF1", "NKG7", "GNLY", "FCGR3A", "CD3D", "CD3G", "CD3E")
DotPlot(object = combined_seurat, features = nk_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")

MajorTypes_panel=c("CD3D",# T cells
                   "MS4A1", # B cells
                   "IGKC", # Plasmablasts
                   "EPCAM", # Epithelial
                   "MKI67", # Cycling
                   "PECAM1", # Endothelial
                   "DCN", # Fibroblasts
                   "COL1A1", # Fibroblast and endothelial
                   "PDPN", # Fibroblasts
                   "LYZ", # Myeloid
                   "TPSAB1", # Mast
                   "GNLY", # NK and T
                   "XCL1", # NK
                   "FGFBP2") # NK

DotPlot(object = combined_seurat, features = MajorTypes_panel) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")




##################################################################################
#T-cell marker checking and NK-NKT marker visualization for AUCells comparisons
##################################################################################

FeaturePlot(combined_seurat, features = "CD3D") + 
  scale_color_viridis_c(option = "magma")

FeaturePlot(combined_seurat, features = "FGFBP2") + 
  scale_color_viridis_c()

FeaturePlot(combined_seurat, features = "FCGR3A") + 
  scale_color_viridis_c()

FeaturePlot(combined_seurat, features = "NCAM1") + 
  scale_color_viridis_c()

FeaturePlot(combined_seurat, features = "KLRF1") + 
  scale_color_viridis_c()

##########################################
#Specific NK and NKT annotation AUCells
##########################################

T_cell_lineage_markers <- c("CD3D", "CD3E", "CD3G")
NK_cell_lineage_markers <- c("GNLY", "XCL1", "FGFBP2")
general_nk_cell_markers <- c("CD7", "NCAM1", "FCGR3A", "B3GAT1", "KLRF1", "NFKB1", "P2RY8", "AXIN1")

# Combine all gene sets into a list
gene_sets <- list(
  "T_Cell_Lineage" = T_cell_lineage_markers,
  "NK_Cell_Lineage" = NK_cell_lineage_markers,
  "General_NK_Cell" = general_nk_cell_markers
)

# Generate expression matrix 
expr_matrix <- as.matrix(GetAssayData(combined_seurat, layer = "data"))

# Build AUCell rankings
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = TRUE)

# Calculate AUC for the gene sets
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)

# Convert AUC to matrix
auc_matrix <- as.matrix(getAUC(cells_AUC))

# Verify AUC matrix dimensions and content
print(dim(auc_matrix))
print(head(auc_matrix))

# Add AUC scores to the metadata of the Seurat object
T_NK_chemo <- AddMetaData(combined_seurat, metadata = t(auc_matrix))

# Convert AUC to long format 
auc_values <- getAUC(cells_AUC)
auc_values_long <- reshape2::melt(auc_values)


# Plot UMAP colored by AUC scores of each gene set
for (gene_set in names(gene_sets)) {
  umap_plot <- FeaturePlot(T_NK_chemo, features = gene_set) + 
    scale_color_viridis_c() + 
    ggtitle(paste("UMAP of", gene_set, "AUC Scores"))
  
  print(umap_plot)
}


T_NK_chemo$T_Cell_Lineage

FeaturePlot(T_NK_chemo, features = "T_Cell_Lineage") + 
  scale_color_viridis_c(option = "magma") +
  ggtitle(paste("UMAP of", "T_Cell_Lineage", "AUC Scores"))


#################################
#Checking CD45 (immune cell) status
################################

DotPlot(object = combined_seurat, features = "PTPRC", group.by = "seurat_clusters") + scale_colour_gradient2(low = "blue", mid = "white", high = "red")

FeaturePlot(combined_seurat, 
            reduction = "umap", 
            features = c("PTPRC"),
            min.cutoff = 'q10', 
            label = TRUE)

############################
#Define immune classification for proportional analyses

cluster_immune_classification <- c(
  "immune",    # 0
  "immune",    # 1
  "non-immune",# 2
  "immune",    # 3
  "immune",    # 4
  "immune",    # 5
  "immune",    # 6
  "immune",    # 7
  "immune",    # 8
  "immune",    # 9
  "non-immune",# 10
  "immune",    # 11
  "immune",    # 12
  "immune",    # 13
  "immune",    # 14
  "immune",    # 15
  "non-immune",# 16
  "immune",    # 17
  "immune",    # 18
  "immune",    # 19
  "immune",    # 20
  "immune",    # 21
  "non-immune",    # 22
  "immune",    # 23
  "immune",    # 24
  "non-immune",    # 25
  "immune",    # 26
  "immune",    # 27
  "immune",    # 28
  "immune",    # 29
  "immune"     # 30
)


combined_seurat$immune_classification <- cluster_immune_classification[combined_seurat$seurat_clusters ]

table(combined_seurat$immune_classification)

DimPlot(combined_seurat, group.by = "immune_classification", reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE)


# Subset 
Pure_immune_seurat <- subset(combined_seurat, subset = immune_classification == 'immune')

# Save 
#saveRDS(Pure_immune_seurat, "Pure_immune_seurat.rds")

###########################################


# Load the data
Pure_immune_seurat <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/Pure_immune_seurat.rds")

#Set back to SCT
#DefaultAssay(combined_seurat) <- "SCT"

##########################################


#######################################
#Proportion analysis
#######################################

Pure_immune_seurat$tissue_tumor <- paste(Pure_immune_seurat$tissue, Pure_immune_seurat$tumorType ,sep = "_")

# Calculate proportions of each cluster within each sample
prop_data <- Pure_immune_seurat@meta.data %>%
  group_by(sample_ID, cluster = Idents(Pure_immune_seurat)) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample_ID) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total) 

prop_data_chemo <- prop_data %>%
  left_join(Pure_immune_seurat@meta.data %>% select(sample_ID, chemo.y), by = "sample_ID")

##############
prop_data_tissue <- prop_data %>%
  left_join(Pure_immune_seurat@meta.data %>% select(sample_ID, tissue), by = "sample_ID")

prop_data_tumor <- prop_data %>%
  left_join(Pure_immune_seurat@meta.data %>% select(sample_ID, tissue_tumor), by = "sample_ID")


####################
###Visualization proportion
##Fig1
sample_order <- paste0("s", 1:25)  # Creates a vector "S1" to "S25"


ggplot(prop_data, aes(x = factor(sample_ID, levels = sample_order), y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Proportion of All Clusters Across All Samples", x = "Sample ID", y = "Proportion")




##Fig2

plot_data_by_category <- function(data, category, category_label) {
  sample_means <- data %>%
    group_by(sample_ID, cluster, !!sym(category)) %>%
    summarise(sample_mean_proportion = mean(proportion), .groups = 'drop')
  
  group_stats <- sample_means %>%
    group_by(cluster, !!sym(category)) %>%
    summarise(
      mean_proportion = mean(sample_mean_proportion),
      se = sd(sample_mean_proportion) / sqrt(n()),  # Standard error
      .groups = 'drop'
    )
  
  dodge_width <- 0.8  
  
  # Plot 
  plot <- ggplot(group_stats, aes(x = cluster, y = mean_proportion, fill = !!sym(category))) +
    geom_col(position = position_dodge(width = dodge_width), color = "black") +
    geom_errorbar(aes(ymin = mean_proportion - se, ymax = mean_proportion + se),
                  width = 0.25 * dodge_width, position = position_dodge(width = dodge_width), color = "black", size = 1) +
    geom_point(data = sample_means, aes(x = cluster, y = sample_mean_proportion),
               position = position_jitterdodge(jitter.width = 0.1, dodge.width = dodge_width), size = 1, color = "black", show.legend = FALSE) + 
    scale_fill_brewer(palette = "Set1", name = category_label) +
    labs(title = paste("Mean Cluster Proportions by", category_label), x = "Cluster", y = "Mean Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plot)
}


plot_chemo <- plot_data_by_category(prop_data_chemo, "chemo.y", "Chemostatus")

filtered_data <- prop_data_tumor %>% 
  filter(grepl("Omentum|Eosophagus|Pancreas", tissue_tumor))

data_omentum <- filtered_data %>% filter(grepl("Omentum", tissue_tumor))
data_esophagus <- filtered_data %>% filter(grepl("Eosophagus", tissue_tumor))
data_pancreas <- filtered_data %>% filter(grepl("Pancreas", tissue_tumor))

plot_omentum <- plot_data_by_category(data_omentum, "tissue_tumor", "Omentum") +
  scale_fill_brewer(palette = "Set1") # Adjusted palette

plot_esophagus <- plot_data_by_category(data_esophagus, "tissue_tumor", "Esophagus") +
  scale_fill_brewer(palette = "Set2") # Adjusted palette

plot_pancreas <- plot_data_by_category(data_pancreas, "tissue_tumor", "Pancreas") +
  scale_fill_brewer(palette = "Set3") # Adjusted palette

combined_plot <- plot_omentum + plot_esophagus + plot_pancreas +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

print(combined_plot)


####################
######
#Proportion analysis
######

# Function to clean data, aggregate proportions, and perform pairwise Wilcoxon test
perform_pairwise_analysis <- function(data, category) {
  data_unique <- data %>%
    distinct(sample_ID, cluster, !!sym(category), .keep_all = TRUE)
  
  data_for_analysis <- data_unique %>%
    group_by(sample_ID, cluster, !!sym(category)) %>%
    summarise(mean_proportion = mean(proportion, na.rm = TRUE), .groups = 'drop')
  
  results_list <- list()
  clusters <- unique(data_for_analysis$cluster)
  
  for (cluster in clusters) {
    data_subset <- data_for_analysis %>% filter(cluster == !!cluster)
    
    if (length(unique(data_subset[[category]])) > 1) {
      pairwise_results <- pairwise_wilcox_test(data_subset, formula = as.formula(paste("mean_proportion ~", category)), p.adjust.method = "BH")
      results_list[[as.character(cluster)]] <- pairwise_results
    } else {
      results_list[[as.character(cluster)]] <- data.frame(
        group1 = NA,
        group2 = NA,
        p = NA,
        p.adj = NA,
        p.format = NA,
        p.signif = NA,
        message = "Insufficient data for comparison in one or more groups"
      )
    }
  }
  
  return(results_list)
}



# Perform analysis for each category
results_chemo <- perform_pairwise_analysis(prop_data_chemo, "chemo.y")
results_tissue <- perform_pairwise_analysis(prop_data_tissue, "tissue")
results_tumor <- perform_pairwise_analysis(prop_data_tumor, "tissue_tumor")

print("Chemostatus Results:")
print(results_chemo)

print("Tissue Type Results:")
print(results_tissue)

print("Tumor Microenvironment Results:")
print(results_tumor)

# Convert results lists to data frames
convert_results_to_df <- function(results_list) {
  # Flatten the list into a data frame
  results_df <- bind_rows(results_list, .id = "cluster")
  return(results_df)
}

# Convert each results list to a data frame
results_chemo_df <- convert_results_to_df(results_chemo)
results_tissue_df <- convert_results_to_df(results_tissue)
results_tumor_df <- convert_results_to_df(results_tumor)

results_chemo_df <- results_chemo_df %>%
  mutate(category_type = "chemo.y")

results_tissue_df <- results_tissue_df %>%
  mutate(category_type = "tissue")

results_tumor_df <- results_tumor_df %>%
  mutate(category_type = "tissue_tumor")

# Combine all results
final_results_df <- bind_rows(results_chemo_df, results_tissue_df, results_tumor_df)

write.csv(final_results_df, file = "/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/Figure_thesis/All_stat_prop_finalized.csv", row.names = FALSE)


