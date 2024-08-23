
##############################################
##########TIDY:HEATMAP#########################
###############################################

### Provide input based on the desired target for visualization

##########################################
##############Chemo################
##########################################
#Prep heatmap data
final_heatmap_data <- final_plot_data

final_heatmap_data <- final_heatmap_data %>%
  mutate(
    NKannotation = str_extract(active_ident, "^[^_]*"),  
#    tumorType = str_extract(active_ident, "(?<=_)[^_]*(?=_)"),  
    chemo = str_extract(active_ident, "[^_]*$")  
  )

final_heatmap_data$variable <- sub("REACTOME_", "", final_heatmap_data$variable)
final_heatmap_data |>
  group_by(NKannotation) |>
  heatmap(
    .column = active_ident,
    .row = variable,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    .value = Mean_AUC,
    scale = "row",
    show_heatmap_legend = TRUE
  ) |> 
  annotation_tile(chemo, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(NKannotation, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) 


############################
#Regulon

final_heatmap_data2 <- final_plot_data2
final_heatmap_data2 <- final_heatmap_data2 %>%
  mutate(
    NKannotation = str_extract(active_ident, "^[^_]*"),  
#    tumorType = str_extract(active_ident, "(?<=_)[^_]*(?=_)"),  
    chemo = str_extract(active_ident, "[^_]*$")  
  )


final_heatmap_data2 <- final_heatmap_data2[order(final_heatmap_data2$NKannotation,final_heatmap_data2$chemo ),] 
  

final_heatmap_data2 |>
  group_by(NKannotation) |>
  heatmap(
    .column = active_ident,
    .row = variable,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    .value = Mean_AUC,
    scale = "row",
    show_heatmap_legend = TRUE
  ) |> 
  annotation_tile(chemo, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(NKannotation, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) 

#####################################################
####Tumor#############
#####################################################


final_heatmap_data <- final_plot_data

final_heatmap_data <- final_heatmap_data %>%
  mutate(
    NKannotation = str_extract(active_ident, "^[^_]*"), 
        tissue = str_extract(active_ident, "(?<=_)[^_]*(?=_)"),  
    tumorType = str_extract(active_ident, "[^_]*$")  
  )

final_heatmap_data$variable <- sub("REACTOME_", "", final_heatmap_data$variable)

final_heatmap_data |>
  group_by(active_ident) |>
  heatmap(
    .column = active_ident,
    .row = variable,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    .value = Mean_AUC,
    scale = "row",
    show_heatmap_legend = TRUE
  ) |> 
    annotation_tile(tumorType, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(tissue, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(NKannotation, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) 


############################
#Regulon

final_heatmap_data2 <- final_plot_data2

final_heatmap_data2 <- final_heatmap_data2 %>%
  mutate(
    NKannotation = str_extract(active_ident, "^[^_]*"),
    tissue = str_extract(active_ident, "(?<=_)[^_]*(?=_)"), 
    tumorType = str_extract(active_ident, "[^_]*$")  
  )

final_heatmap_data2$variable <- sub("GO:BP_", "", final_heatmap_data2$variable)


final_heatmap_data2 |>
  group_by(active_ident) |>
  heatmap(
    .column = active_ident,
    .row = variable,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    .value = Mean_AUC,
    scale = "row",
    show_heatmap_legend = TRUE
  ) |> 
   annotation_tile(tumorType, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(tissue, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(NKannotation, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) 









#####################################################
####TumorChemo#############
#####################################################

final_heatmap_data <- final_plot_data

final_heatmap_data <- final_heatmap_data %>%
  mutate(
    split_ident = str_split(active_ident, "_", simplify = TRUE),
    NKannotation = split_ident[,1],  
    tissue = split_ident[,2], 
    tumorType = split_ident[,3],  
    chemo = split_ident[,4] 
  ) %>%
  select(-split_ident)  


final_heatmap_data$variable <- sub("REACTOME_", "", final_heatmap_data$variable)


final_heatmap_data |>
  group_by(active_ident) |>
  heatmap(
    .column = active_ident,
    .row = variable,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    .value = Mean_AUC,
    scale = "row",
    show_heatmap_legend = TRUE
  ) |> 
  annotation_tile(tumorType, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(tissue, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(chemo, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |>
  annotation_tile(NKannotation, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) 




##################################
#Regulon
##################################

final_heatmap_data2 <- final_plot_data2

final_heatmap_data2 <- final_heatmap_data2 %>%
  mutate(
    split_ident = str_split(active_ident, "_", simplify = TRUE),
    NKannotation = split_ident[,1],  
    tissue = split_ident[,2], 
    tumorType = split_ident[,3],  
    chemo = split_ident[,4]  
  ) %>%
  select(-split_ident) 


final_heatmap_data2 |>
  group_by(active_ident) |>
  heatmap(
    .column = active_ident,
    .row = variable,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    .value = Mean_AUC,
    scale = "row",
    show_heatmap_legend = TRUE
  ) |> 
  annotation_tile(tumorType, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |>
  annotation_tile(tissue, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(chemo, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) |> 
  annotation_tile(NKannotation, size = unit(0.3, "cm"),  annotation_name_gp= gpar(fontsize = 8)) 

###########################################################################

###########################################################################
#AUCells UMAP comparisons
#############################################################################

###ONE OF THE EXAMPLE 
#The rest were manually changing input code

AUC_pathway_seurat_UMAP_tissue <- AUC_regulon_seurat
AUC_pathway_seurat_UMAP_tissue$NK_tissue_tumor <- paste(AUC_pathway_seurat_UMAP_tissue$tissue, AUC_pathway_seurat_UMAP_tissue$tumorType ,sep = "_")

Idents(AUC_pathway_seurat_UMAP_tissue) <- "NK_tissue_tumor"
print(table(Idents(AUC_pathway_seurat_UMAP_tissue)))

AUC_pathway_seurat_targettissue <- subset(AUC_pathway_seurat_UMAP_tissue, idents = c("Eosophagus_Norm"))

AUC_pathway_seurat_others <- subset(AUC_pathway_seurat_UMAP_tissue, idents = c("Eosophagus_Norm"), invert = TRUE)

AUC_pathway_seurat_UMAP_tissue@meta.data$LAMB3_TARGET_GENES

UMAP_target <- FeaturePlot(AUC_pathway_seurat_targettissue, features = "LAMB3_TARGET_GENES") + 
  scale_color_viridis_c()

UMAP_others <- FeaturePlot(AUC_pathway_seurat_others, features = "LAMB3_TARGET_GENES") + 
  scale_color_viridis_c()

UMAP_tissue_combined <- UMAP_target + UMAP_others

#############################################################################
