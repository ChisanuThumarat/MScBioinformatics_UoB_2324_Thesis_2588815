
#Set back to SCT
DefaultAssay(T_NK_cell_subset) <- "SCT"

##############################
#First processing after subset
##############################

DimPlot(T_NK_cell_subset, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(T_NK_cell_subset, reduction = "umap", group.by = "seurat_clusters",label = TRUE, pt.size = 0.5) 

######################################
T_NK_cell_subset <- FindVariableFeatures(T_NK_cell_subset, selection.method = "vst", nfeatures = 2000, assay="RNA")

#Mitochondrial percentage for SCT regression 
T_NK_cell_subset <- PercentageFeatureSet(T_NK_cell_subset, pattern = "^MT-", col.name = "percent.mt")

#SCTransform normalization
T_NK_cell_subset <- SCTransform(T_NK_cell_subset, vars.to.regress = "percent.mt", verbose = FALSE)

#Rerun Harmony for batch correction
T_NK_cell_subset <- RunPCA(T_NK_cell_subset, assay.use = "SCT",  verbose = FALSE)
T_NK_cell_subset <- RunHarmony(T_NK_cell_subset, assay.use = "SCT",  "sample_ID")
Eplot<-ElbowPlot(T_NK_cell_subset)
T_NK_cell_subset  <- RunUMAP(T_NK_cell_subset, assay = "SCT",  dims = 1:20) # Create new UMAP


#UMAP
DimPlot(T_NK_cell_subset, reduction = "umap", group.by = "sample_ID") + ggtitle("UMAP Plot Combined Seurat")
#Harmony
DimPlot(T_NK_cell_subset, reduction = "harmony", group.by = "sample_ID") + ggtitle("Harmony Plot Combined Seurat")

# Find Neighbors
T_NK_cell_subset <- FindNeighbors(T_NK_cell_subset, assay = "SCT", reduction = "harmony", dims = 1:20)
# Find Clusters
T_NK_cell_subset <- FindClusters(T_NK_cell_subset, resolution = 0.5, graph.name = "SCT_snn")
DimPlot(T_NK_cell_subset, reduction = "umap", label = TRUE)


#saveRDS(T_NK_cell_subset, "T_NK_cell_Ready_to_cluster.rds")

#T_NK_cell_subset <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/T_NK_cell_Ready_to_cluster.rds")


####################################
#Check and try to isolate NK and NKT
####################################

#Dotplot
nk_markers <- c("NCAM1", "KLRF1", "NKG7", "GNLY", "FCGR3A", "CD3D", "CD3G", "CD3E")
DotPlot(object = T_NK_cell_subset, features = nk_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")



NK_cell_CD56 <- FeaturePlot(T_NK_cell_subset, 
                            reduction = "umap", 
                            features = c("NCAM1"),
                            min.cutoff = 'q10', 
                            label = TRUE)

NK_cell_NKG7 <- FeaturePlot(T_NK_cell_subset, 
                            reduction = "umap", 
                            features = c("NKG7"),
                            min.cutoff = 'q10', 
                            label = TRUE)

NK_cell_GNLY <- FeaturePlot(T_NK_cell_subset, 
                            reduction = "umap", 
                            features = c("GNLY"),
                            min.cutoff = 'q10', 
                            label = TRUE)


NK_cell_KLRF1 <- FeaturePlot(T_NK_cell_subset, 
                             reduction = "umap", 
                             features = c("KLRF1"),
                             min.cutoff = 'q10', 
                             label = TRUE)

NK_cell_FCGR3A <- FeaturePlot(T_NK_cell_subset, 
                              reduction = "umap", 
                              features = c("FCGR3A"),
                              min.cutoff = 'q10', 
                              label = TRUE)

NK_cell_CD3D <- FeaturePlot(T_NK_cell_subset, 
                            reduction = "umap", 
                            features = c("CD3D"),
                            min.cutoff = 'q10', 
                            label = TRUE)

NK_cell_CD3G <- FeaturePlot(T_NK_cell_subset, 
                            reduction = "umap", 
                            features = c("CD3G"),
                            min.cutoff = 'q10', 
                            label = TRUE)

NK_cell_CD3E <- FeaturePlot(T_NK_cell_subset, 
                            reduction = "umap", 
                            features = c("CD3E"),
                            min.cutoff = 'q10', 
                            label = TRUE)
NK_cell_CD56
NK_cell_KLRF1
NK_cell_FCGR3A
NK_cell_CD3D 
NK_cell_CD3E 
NK_cell_CD3G 

new.cluster.ids <- c(
  0,  # 0
  1,  # 1
  'NK',  # 2
  3,  # 3
  4,  # 4
  5,  # 5
  'NK',  # 6
  7,  # 7
  'NK',  # 8
  'NK',  # 9
  10,  # 10
  11,  # 11
  12,  # 12
  'NK',  # 13
  'NK'  # 14
)
names(new.cluster.ids) <- levels(T_NK_cell_subset)
T_NK_cell_subset <- RenameIdents(T_NK_cell_subset, new.cluster.ids)
DimPlot(T_NK_cell_subset, reduction = "umap", label = TRUE, pt.size = 0.5) 


#Subset
T_NK_2nd_subset <- subset(T_NK_cell_subset, idents = "NK")
DimPlot(T_NK_2nd_subset, reduction = "umap", label = TRUE, pt.size = 0.5) 


# Save after the second round 
#saveRDS(T_NK_2nd_subset, "T_NK_2nd_subset.rds")
###############################################################################


# Load the data
#T_NK_2nd_subset <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/T_NK_2nd_subset.rds")

#Set back to SCT
#DefaultAssay(T_NK_2nd_subset) <- "SCT"

################################################################################

##############################
#Second processing after subset
##############################

T_NK_2nd_subset <- FindVariableFeatures(T_NK_2nd_subset, selection.method = "vst", nfeatures = 2000, assay="RNA")
T_NK_2nd_subset <- PercentageFeatureSet(T_NK_2nd_subset, pattern = "^MT-", col.name = "percent.mt")
T_NK_2nd_subset <- SCTransform(T_NK_2nd_subset, vars.to.regress = "percent.mt", verbose = FALSE)
T_NK_2nd_subset <- RunPCA(T_NK_2nd_subset, assay.use = "SCT",  verbose = FALSE)
T_NK_2nd_subset <- RunHarmony(T_NK_2nd_subset, assay.use = "SCT",  "sample_ID")
Eplot<-ElbowPlot(T_NK_2nd_subset)
T_NK_2nd_subset  <- RunUMAP(T_NK_2nd_subset, assay = "SCT", dims = 1:20)


#UMAP
DimPlot(T_NK_2nd_subset, reduction = "umap", group.by = "sample_ID") + ggtitle("UMAP Plot Combined Seurat")
#Harmony
DimPlot(T_NK_2nd_subset, reduction = "harmony", group.by = "sample_ID") + ggtitle("Harmony Plot Combined Seurat")

# Find Neighbors
T_NK_2nd_subset <- FindNeighbors(T_NK_2nd_subset, reduction = "harmony", dims = 1:20)
# Find Clusters
T_NK_2nd_subset <- FindClusters(T_NK_2nd_subset, resolution = 0.5)
DimPlot(T_NK_2nd_subset, reduction = "umap", label = TRUE)


#########################
#Isolation step
########################
#Dotplot
nk_markers <- c("NCAM1", "KLRF1", "NKG7", "GNLY", "FCGR3A", "CD3D", "CD3G", "CD3E")
DotPlot(object = T_NK_2nd_subset, features = nk_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")

NK_diff=c("CD7", "NCAM1","B3GAT1","KLRD1", "FCGR3A", "NKG7", "GNLY","IL7R","SELL","KLRC1","KLRC2","PDCD1","HAVCR2","CD69","NCR1","KLRK1","CD44","CD226","XCL1","XCL2","CD96","MKI67")
DotPlot(object = T_NK_2nd_subset, features = NK_diff) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")


FeaturePlot(T_NK_2nd_subset, 
            reduction = "umap", 
            features = c("NCAM1"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_2nd_subset, 
            reduction = "umap", 
            features = c("NKG7"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_2nd_subset, 
            reduction = "umap", 
            features = c("GNLY"),
            min.cutoff = 'q10', 
            label = TRUE)


FeaturePlot(T_NK_2nd_subset, 
            reduction = "umap", 
            features = c("KLRF1"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_2nd_subset, 
            reduction = "umap", 
            features = c("FCGR3A"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_2nd_subset, 
            reduction = "umap", 
            features = c("CD3D"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_2nd_subset, 
            reduction = "umap", 
            features = c("CD3G"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_2nd_subset, 
            reduction = "umap", 
            features = c("CD3E"),
            min.cutoff = 'q10', 
            label = TRUE)

new.cluster.ids <- c(
  'NK',  # 0
  1,  # 1
  'NK',  # 2
  'NK',  # 3
  'NK',  # 4
  'NK',  # 5
  'NK',  # 6
  'NK',  # 7
  'NK',  # 8
  'NK',  # 9
  10 , # 10
  'NK'
)
names(new.cluster.ids) <- levels(T_NK_2nd_subset)
T_NK_2nd_subset <- RenameIdents(T_NK_2nd_subset, new.cluster.ids)
DimPlot(T_NK_2nd_subset, reduction = "umap", label = TRUE, pt.size = 0.5) 

#Dotplot
nk_markers <- c("NCAM1", "KLRF1", "NKG7", "GNLY", "FCGR3A", "CD3D", "CD3G", "CD3E")
DotPlot(object = T_NK_2nd_subset, features = nk_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# Subset the Seurat object for T_NK_cell identity
T_NK_3rd_subset <- subset(T_NK_2nd_subset, idents = "NK")

##Save 3rd round
saveRDS(T_NK_3rd_subset, "T_NK_3rd_subset.rds")
######################################################

# Load the data
#T_NK_3rd_subset <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/T_NK_3rd_subset.rds")

#Set back to SCT
#DefaultAssay(T_NK_3rd_subset) <- "SCT"

#######################################################


##############################
#Third processing after subset
##############################
T_NK_3rd_subset <- FindVariableFeatures(T_NK_3rd_subset, selection.method = "vst", nfeatures = 2000, assay="RNA")
T_NK_3rd_subset <- PercentageFeatureSet(T_NK_3rd_subset, pattern = "^MT-", col.name = "percent.mt")
T_NK_3rd_subset <- SCTransform(T_NK_3rd_subset, vars.to.regress = "percent.mt", verbose = FALSE)
T_NK_3rd_subset <- RunPCA(T_NK_3rd_subset, assay.use = "SCT",  verbose = FALSE)
T_NK_3rd_subset <- RunHarmony(T_NK_3rd_subset, assay.use = "SCT",  "sample_ID")
Eplot<-ElbowPlot(T_NK_3rd_subset)
T_NK_3rd_subset  <- RunUMAP(T_NK_3rd_subset, assay = "SCT", dims = 1:20) # Create new UMAP
DimPlot(T_NK_3rd_subset, reduction = "harmony", group.by = "sample_ID") + ggtitle("Harmony Plot Combined Seurat")

# Find Neighbors
T_NK_3rd_subset <- FindNeighbors(T_NK_3rd_subset, reduction = "harmony", dims = 1:20)
# Find Clusters
T_NK_3rd_subset <- FindClusters(T_NK_3rd_subset, resolution = 0.5)
DimPlot(T_NK_3rd_subset, reduction = "umap", label = TRUE)


#########################
#Isolation step
########################

#Dotplot
nk_markers <- c("NCAM1", "KLRF1", "NKG7", "GNLY", "FCGR3A", "CD3D", "CD3G", "CD3E")
DotPlot(object = T_NK_3rd_subset, features = nk_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")

NK_diff=c("CD7", "NCAM1","B3GAT1","KLRD1", "FCGR3A", "NKG7", "GNLY","IL7R","SELL","KLRC1","KLRC2","PDCD1","HAVCR2","CD69","NCR1","KLRK1","CD44","CD226","XCL1","XCL2","CD96","MKI67")
DotPlot(object = T_NK_3rd_subset, features = NK_diff) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")


FeaturePlot(T_NK_3rd_subset, 
            reduction = "umap", 
            features = c("NCAM1"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_3rd_subset, 
            reduction = "umap", 
            features = c("NKG7"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_3rd_subset, 
            reduction = "umap", 
            features = c("GNLY"),
            min.cutoff = 'q10', 
            label = TRUE)


FeaturePlot(T_NK_3rd_subset, 
            reduction = "umap", 
            features = c("KLRF1"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_3rd_subset, 
            reduction = "umap", 
            features = c("FCGR3A"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_3rd_subset, 
            reduction = "umap", 
            features = c("CD3D"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_3rd_subset, 
            reduction = "umap", 
            features = c("CD3G"),
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(T_NK_3rd_subset, 
            reduction = "umap", 
            features = c("CD3E"),
            min.cutoff = 'q10', 
            label = TRUE)

RidgePlot(object = T_NK_3rd_subset, features = 'CD3D')
RidgePlot(object = T_NK_3rd_subset, features = 'NCAM1')
RidgePlot(object = T_NK_3rd_subset, features = 'KLRF1')
RidgePlot(object = T_NK_3rd_subset, features = 'NKG7')
RidgePlot(object = T_NK_3rd_subset, features = 'FCGR3A')

##################################################
# Differential expression analysis for clusters
##################################################

DefaultAssay(T_NK_3rd_subset) <- "RNA"
NormalizeData(T_NK_3rd_subset,normalization.method = "LogNormalize", scale.factor = 10000)
FindVariableFeatures(T_NK_3rd_subset, selection.method = "vst", nfeatures = 2000, assay="RNA")
ScaleData(T_NK_3rd_subset)
Nk.markers <- FindAllMarkers(JoinLayers(T_NK_3rd_subset),assay="RNA",  min.pct = 0.25, logfc.threshold = 0.25)
DefaultAssay(T_NK_3rd_subset) <- "SCT"

# To find the top 10 markers
top_markers <- Nk.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(T_NK_3rd_subset, features = top_markers$gene)



#####################################
#Finalized grouping checking 
####################################

new.cluster.ids <- c(
  'NKg1',  # 0
  'NKg2',  # 1
  'NKg3',  # 2
  'NKg1',  # 3
  'NKg4',  # 4
  'NKg5',  # 5
  'NKg6',  # 6
  'Unknown',  # 7
  'Unknown',  # 8
  'NKg7',  # 9
  'Unknown',  # 10
  'NKg8'  # 11
)

names(new.cluster.ids) <- levels(T_NK_3rd_subset)
T_NK_3rd_subset <- RenameIdents(T_NK_3rd_subset, new.cluster.ids)
DimPlot(T_NK_3rd_subset, reduction = "umap", label = TRUE, pt.size = 0.5) 

##Save finalized version 
saveRDS(T_NK_3rd_subset, "T_NK_finalized_subset.rds")


######################################################

# Load the data
T_NK_3rd_subset <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/T_NK_finalized_subset.rds")

##Re-labelled after checking NK groups 

###############################
#Finalized NK data based on literature and identity
##############################
# Mapping from old labels to new labels
label.mapping <- c(
  'NKg1' = 'NKcyt',
  'NKg2' = 'NKexh',
  'NKg3' = 'NKimm',
  'NKg4' = 'NKeff1',
  'NKg5' = 'NKeff2',
  'NKg6' = 'NKmem',
  'NKg7' = 'NKres',
  'NKg8' = 'NKact',
  'Unknown' = 'Unknown'
)


current.cluster.ids <- factor(Idents(T_NK_3rd_subset))
new.cluster.ids <- plyr::mapvalues(current.cluster.ids, from = names(label.mapping), to = label.mapping)
Idents(T_NK_3rd_subset) <- new.cluster.ids
print(Idents(T_NK_3rd_subset))

saveRDS(T_NK_3rd_subset, "T_NK_finalized_relabelled_subset.rds")

####################################

# Load the data
T_NK_3rd_subset <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/T_NK_finalized_relabelled_subset.rds")


####################################
#Functional NK markers testing
####################################
#Marker to test
cytotoxicity_markers <- c("GNLY", "GZMA", "GIMAP7", "PSMB10", "TBC1D10C", "FGFBP2", "ACTB", "CSK")
activation_markers <- c("CD69", "NCR1", "KLRK1", "CD226", "CD96", "FOS", "IFNG", "ATF3", "PELI2", "IGF1R", "IL12RB2", "IER3", "EGR1", "AREG")
adhesion_and_migration_markers <- c("CD44", "SELL", "CXCR6", "CHPT1", "TIAM1")
proliferation_markers <- c("MKI67", "PVT1", "PBX4", "IRS2", "IL7R")
inhibitory_markers <- c("KLRD1", "KLRC1", "KLRC2", "PDCD1", "HAVCR2", "KIR2DL4", "BACH2", "NFKBID", "IKZF2")
chemokines_and_receptors_markers <- c("XCL1", "XCL2", "CCL3", "CCL4", "CCL4L2", "CCL3L1")
stress_and_exhaustion_markers <- c("HSPA1A", "HSP90AA1", "DNAJB1", "HSPA1B", "HSPH1", "HSPD1", "FKBP4", "HSPA6", "HSPB1", "ZFAND2A")
general_nk_cell_markers <- c("CD7", "NCAM1", "FCGR3A", "B3GAT1", "KLRF1", "NFKB1", "P2RY8", "PHLDB2", "BNC2", "AXIN1", "INPP5A")
miscellaneous_other_functional_roles <- c("CRTAM", "LINC02446", "RILPL2", "MAML3", "B3GNT7", "SSBP2", "ABTB2", "FTH1", "C1orf21", "TGFBR3", "KLF3", "PEX14", "OASL", "WDR47", "VPS13D")
features_to_test <- list(
  "Cytotoxicity" = cytotoxicity_markers,
  "Adhesion_and_Migration" = adhesion_and_migration_markers
)

features_to_test2 <- list(
  "Proliferation" = proliferation_markers,
  "Inhibitory" = inhibitory_markers,
  "Chemokines_and_Receptors" = chemokines_and_receptors_markers
)

features_to_test3 <- list(
  "Stress_and_Exhaustion" = stress_and_exhaustion_markers,
  "General_NK_Cell" = general_nk_cell_markers
  
)


features_to_test4 <- list(
  "Activation" = activation_markers
)

features_to_test5 <- list(
  "Miscellaneous_Other_Functional_Roles" = miscellaneous_other_functional_roles
)


DotPlot(object = T_NK_3rd_subset, features = features_to_test) + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") 

DotPlot(object = T_NK_3rd_subset, features = features_to_test2) + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") 

DotPlot(object = T_NK_3rd_subset, features = features_to_test3) + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") 

DotPlot(object = T_NK_3rd_subset, features = features_to_test4) + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") 

DotPlot(object = T_NK_3rd_subset, features = features_to_test5) + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") 



#Dotplot

DotPlot(object = T_NK_3rd_subset, features = cytotoxicity_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red") + RotatedAxis()+coord_flip()

DotPlot(object = T_NK_3rd_subset, features = activation_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+ RotatedAxis()+coord_flip()

DotPlot(object = T_NK_3rd_subset, features = adhesion_and_migration_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+ RotatedAxis()+coord_flip()

DotPlot(object = T_NK_3rd_subset, features = proliferation_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+ RotatedAxis()+coord_flip()

DotPlot(object = T_NK_3rd_subset, features = inhibitory_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+ RotatedAxis()+coord_flip()

DotPlot(object = T_NK_3rd_subset, features = chemokines_and_receptors_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+ RotatedAxis()+coord_flip()

DotPlot(object = T_NK_3rd_subset, features = stress_and_exhaustion_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+ RotatedAxis()+coord_flip()

DotPlot(object = T_NK_3rd_subset, features = general_nk_cell_markers) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+ RotatedAxis()+coord_flip()

DotPlot(object = T_NK_3rd_subset, features = miscellaneous_other_functional_roles) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+ RotatedAxis()+coord_flip()



#########################################
#AUCells-based on fucntional markers
#########################################
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


expr_matrix <- as.matrix(GetAssayData(T_NK_3rd_subset, layer = "data"))

if (is.null(expr_matrix) || ncol(expr_matrix) == 0 || nrow(expr_matrix) == 0) {
  stop("The expression matrix is empty. Please check the specified layer and data extraction.")
}


###AUCells ranking stage
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = TRUE)

cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)

auc_matrix <- as.matrix(getAUC(cells_AUC))

# AUC assay to store information
T_NK_3rd_subset <- AddMetaData(T_NK_3rd_subset, metadata = t(auc_matrix))


#UMAP AUCells
for (gene_set in names(gene_sets)) {
  umap_plot <- FeaturePlot(T_NK_3rd_subset, features = gene_set) + 
    scale_color_viridis_c() + 
    ggtitle(paste("UMAP of", gene_set, "AUC Scores"))
  
  print(umap_plot)
}



#######################################
#Proportion analysis
#######################################

T_NK_3rd_subset$NK_tissue_tumor <- paste(T_NK_3rd_subset$tissue, T_NK_3rd_subset$tumorType ,sep = "_")

# Calculate proportions of each cluster within each sample
prop_data <- T_NK_3rd_subset@meta.data %>%
  group_by(sample_ID, cluster = Idents(T_NK_3rd_subset)) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample_ID) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total)

### Prep the data by conditions
prop_data_tissue <- prop_data %>%
  left_join(T_NK_3rd_subset@meta.data %>% select(sample_ID, tissue), by = "sample_ID")

prop_data_chemo <- prop_data %>%
  left_join(T_NK_3rd_subset@meta.data %>% select(sample_ID, chemo.y), by = "sample_ID")

##############
prop_data_tumor <- prop_data %>%
  left_join(T_NK_3rd_subset@meta.data %>% select(sample_ID, NK_tissue_tumor), by = "sample_ID")


####################

#Visualization

##Fig1
sample_order <- paste0("s", 1:25) 

ggplot(prop_data, aes(x = factor(sample_ID, levels = sample_order), y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Proportion of All NK functional states Across All Samples", x = "Sample ID", y = "Proportion")

##Fig2

plot_data_by_category <- function(data, category, category_label) {
  # Calculate sample-level mean proportions
  sample_means <- data %>%
    group_by(sample_ID, cluster, !!sym(category)) %>%
    summarise(sample_mean_proportion = mean(proportion), .groups = 'drop')
  
  # Calculate group-level statistics
  group_stats <- sample_means %>%
    group_by(cluster, !!sym(category)) %>%
    summarise(
      mean_proportion = mean(sample_mean_proportion),
      se = sd(sample_mean_proportion) / sqrt(n()),  # Standard error
      .groups = 'drop'
    )
  
  
  dodge_width <- 0.8  
  
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


#######################


filtered_data <- prop_data_tumor %>% 
  filter(grepl("Omentum|Eosophagus|Pancreas", NK_tissue_tumor))

# Then, split the data by the main categories to apply the function separately
data_omentum <- filtered_data %>% filter(grepl("Omentum", NK_tissue_tumor))
data_esophagus <- filtered_data %>% filter(grepl("Eosophagus", NK_tissue_tumor))
data_pancreas <- filtered_data %>% filter(grepl("Pancreas", NK_tissue_tumor))

# Apply the plotting function to each subset with unique color palettes
plot_omentum <- plot_data_by_category(data_omentum, "NK_tissue_tumor", "Omentum") +
  scale_fill_brewer(palette = "Set1") # Adjusted palette

plot_esophagus <- plot_data_by_category(data_esophagus, "NK_tissue_tumor", "Esophagus") +
  scale_fill_brewer(palette = "Set2") # Adjusted palette

plot_pancreas <- plot_data_by_category(data_pancreas, "NK_tissue_tumor", "Pancreas") +
  scale_fill_brewer(palette = "Set3") # Adjusted palette

# Arrange the plots into a single visual
combined_plot <- plot_omentum + plot_esophagus + plot_pancreas +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Print the combined plot
print(combined_plot)

########
plot_chemo <- plot_data_by_category(prop_data_chemo, "chemo.y", "Chemostatus")



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

##Prepare data for analysis
prop_data <- T_NK_3rd_subset@meta.data %>%
  group_by(sample_ID, cluster = Idents(T_NK_3rd_subset)) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample_ID) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total)

#########################################################
# Perform analysis for each category
results_chemo <- perform_pairwise_analysis(prop_data_chemo, "chemo.y")
results_tissue <- perform_pairwise_analysis(prop_data_tissue, "tissue")
results_tumor <- perform_pairwise_analysis(prop_data_tumor, "NK_tissue_tumor")


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
  mutate(category_type = "tumorType")

# Combine all results into a single data frame
final_results_df <- bind_rows(results_chemo_df, results_tissue_df, results_tumor_df)


################################################


write.csv(final_results_df, file = "/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/Figure_thesis/All_stat_prop_NK_finalized.csv", row.names = FALSE)


