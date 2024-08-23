# Load the data
Seurat_object_list <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/DATA/xxx")
Meta_data <- read.csv("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/DATA/xxx")
#Set seed
set.seed(42) 

############################
#QC processing stored in list
############################

# Initialize a new list to store processed Seurat objects
processed_seurat_objects <- list()

for (i in seq_along(Seurat_object_list)) {
  seurat_object <- Seurat_object_list[[i]]
  
  # Calculate the percentage of mitochondrial genes 
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
  sum_stat_qc_plot  <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
    ggtitle(paste("Sample", i)) 
  ggsave(paste0("sum_stat_preqc", i, ".png"), plot = sum_stat_qc_plot, width = 8, height = 6)
  
  plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")+
    ggtitle(paste("Sample", i)) 
  plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    ggtitle(paste("Sample", i)) 
  feature_scatter_preqc <- plot1 + plot2+
    ggtitle(paste("Sample", i)) 
  ggsave(paste0("feature_scatter_preqc", i, ".png"), plot = feature_scatter_preqc , width = 12, height = 8)
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 300 & nFeature_RNA < 3500 & percent.mt < 10)
  sum_stat_post_qc_plot  <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
    ggtitle(paste("Sample", i)) 
  ggsave(paste0("sum_stat_postqc", i, ".png"), plot = sum_stat_post_qc_plot, width = 8, height = 6)
  
  plot3 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")+
    ggtitle(paste("Sample", i)) 
  plot4 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    ggtitle(paste("Sample", i)) 
  feature_scatter_postqc <- plot3 + plot4+
    ggtitle(paste("Sample", i)) 
  ggsave(paste0("feature_scatter_postqc", i, ".png"), plot = feature_scatter_postqc , width = 12, height = 8)
  
  seurat_object <- SCTransform(seurat_object, vars.to.regress = "percent.mt", verbose = FALSE)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  
  top10 <- head(VariableFeatures(seurat_object), 10)
  
  plot5 <- VariableFeaturePlot(seurat_object)+
    ggtitle(paste("Sample", i)) 
  plot6 <- LabelPoints(plot = plot5, points = top10, repel = TRUE)+
    ggtitle(paste("Sample", i)) 
  vst_predouF <- plot5 + plot6+
    ggtitle(paste("Sample", i)) 
  
  ggsave(paste0("vst_predouF", i, ".png"), plot = vst_predouF , width = 12, height = 8)
  
  processed_seurat_objects[[i]] <- seurat_object
}

##################
#PCA

num_pcs <- 20 
plot_list <- list()

for (i in seq_along(processed_seurat_objects)) {
  processed_seurat_objects[[i]] <- RunPCA(processed_seurat_objects[[i]], features = VariableFeatures(object = processed_seurat_objects[[i]]))
  p <- DimPlot(processed_seurat_objects[[i]], reduction = "pca", group.by = "sample_ID")
  plot_list[[i]] <- p
  ggsave(paste0("PCA_plot_sample_", i, ".png"), plot = p, width = 6, height = 4)
}
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 3)
ggsave("combined_PCA_plots.png", plot = combined_plot, width = 18, height = 12)




########################
#DoubletFinder
########################

num_pcs <- 20 

processed_seurat_objects <- lapply(processed_seurat_objects, function(x) {
  x <- doubletFinder(
    x,
    PCs = 1:num_pcs,  
    pN = 0.25,
    pK = 0.01,
    nExp = ncol(x) * 0.05,
    reuse.pANN = FALSE,
    sct = TRUE 
  )
  
  return(x)
})






