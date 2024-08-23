
T_NK_3rd_subset <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/T_NK_finalized_relabelled_subset.rds")

#####################################
##For chemo and tumor

#Prepare for storing chemostatus
T_NK_chemo  <- T_NK_3rd_subset
T_NK_chemo$NK.chemo <- paste(T_NK_chemo$NK_annotation ,T_NK_chemo$chemo.y, sep = "_")

Idents(T_NK_chemo) <- "NK.chemo"
print(table(Idents(T_NK_chemo)))



T_NK_tumorType  <- T_NK_3rd_subset
T_NK_tumorType$NK.tumor <- paste(T_NK_tumorType$NK_annotation ,T_NK_tumorType$tissue ,T_NK_tumorType$tumorType, sep = "_")

Idents(T_NK_tumorType) <- "NK.tumor"
print(table(Idents(T_NK_tumorType)))


############################

##Case subset for chemotumor

Idents(T_NK_3rd_subset) <- T_NK_3rd_subset@meta.data$tumorType
T_NK_3rd_subset_case_only <- subset(T_NK_3rd_subset, idents = c("GI_MET", "Benign", "PDAC"), invert = TRUE)

all(sapply(T_NK_3rd_subset_case_only@assays, function(x) all(rownames(x) %in% rownames(T_NK_3rd_subset_case_only@assays$RNA))))
common_genes <- Reduce(intersect, lapply(T_NK_3rd_subset_case_only@assays, rownames))

##Intersect for making the matrix retain only the gene of the leftovers
T_NK_3rd_subset_case_only <- T_NK_3rd_subset_case_only[common_genes, ]

gi_met_genes <- rownames(subset(T_NK_3rd_subset, idents = c("GI_MET", "Benign", "PDAC"))@assays$RNA)
other_genes <- rownames(subset(T_NK_3rd_subset, idents =c("GI_MET", "Benign", "PDAC"), invert = TRUE)@assays$RNA)

unique_gi_met_genes <- setdiff(gi_met_genes, other_genes)
unique_other_genes <- setdiff(other_genes, gi_met_genes)

##################################################################

# Verify the gene consistency
consistent_genes_after <- all(sapply(T_NK_3rd_subset_case_only@assays, function(x) all(rownames(x) == common_genes)))
print(paste("Gene consistency across layers after adjustment:", consistent_genes_after))

# Normalize and find variable features
T_NK_3rd_subset_case_only <- NormalizeData(T_NK_3rd_subset_case_only, verbose = FALSE)
T_NK_3rd_subset_case_only <- FindVariableFeatures(T_NK_3rd_subset_case_only, selection.method = "vst", nfeatures = 2000)

# Scale data
T_NK_3rd_subset_case_only <- ScaleData(T_NK_3rd_subset_case_only, features = rownames(T_NK_3rd_subset_case_only), verbose = FALSE)
Idents(T_NK_3rd_subset_case_only) <- T_NK_3rd_subset_case_only@meta.data$NK_annotation
saveRDS(T_NK_3rd_subset_case_only, "T_NK_case_only.rds")

#####################################

T_NK_3rd_subset <- readRDS("/rds/projects/m/mossp-ovarian-cancer/CHARLIE/thesis_project/T_NK_case_only.rds")

#################################

###For now use all cases
T_NK_both  <- T_NK_3rd_subset_case_only
T_NK_both$NK.both <- paste(T_NK_both$NK_annotation ,T_NK_both$tissue, T_NK_both$tumorType,  T_NK_both$chemo.y, sep = "_")

Idents(T_NK_both) <- "NK.both"
print(table(Idents(T_NK_both)))

#saveRDS(T_NK_chemo, "T_NK_chemo_AUCells.rds")
saveRDS(T_NK_tumorType, "T_NK_tumorType_AUCells.rds")
#saveRDS(T_NK_both, "T_NK_both_AUCells.rds")


saveRDS(T_NK_chemo, "T_NK_chemo_AUCells_case_only.rds")
saveRDS(T_NK_tumorType, "T_NK_tumorType_AUCells_case_only.rds")
saveRDS(T_NK_both, "T_NK_both_AUCells_case_only.rds")
