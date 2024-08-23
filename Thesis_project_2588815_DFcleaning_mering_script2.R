
###Run this in case matrix overloaded in R
##Using batch integration and merging is trading off to save progression
#Set seed
set.seed(42) 

##Function 
divide_into_batches <- function(obj_list, batch_size) {
  split(obj_list, ceiling(seq_along(obj_list) / batch_size))
}


batches <- divide_into_batches(cleaned_seurat_objects, 5)

######################################################

## Run Anchors and features
##########################################################
# lists 
features_list <- list()
anchors_list <- list()

for (i in seq_along(batches)) {
  features <- SelectIntegrationFeatures(object.list = batches[[i]], nfeatures = 3000)
  features_list[[i]] <- features
  batches[[i]] <- PrepSCTIntegration(object.list = batches[[i]], anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list = batches[[i]], normalization.method = "SCT", anchor.features = features)
  anchors_list[[i]] <- anchors
  saveRDS(anchors, paste0("anchors_batch_", i, ".rds"))
}


###############
# List : Integrated batch
integrated_batches <- list()

for (i in seq_along(anchors_list)) {
  integrated_data <- IntegrateData(anchorset = anchors_list[[i]], normalization.method = "SCT")
  integrated_batches[[i]] <- integrated_data
  saveRDS(integrated_data, paste0("integrated_batch_", i, ".rds"))
}


##########################
#Continue from Batch file
##########################


# Load file function
load_integrated_batches <- function(batch_indices) {
  integrated_batches <- list()
  for (i in batch_indices) {
    integrated_batch <- readRDS(paste0("integrated_batch_", i, ".rds"))
    integrated_batches[[i]] <- integrated_batch
  }
  return(integrated_batches)
}


#####################################
#Removing DoubletFinder traces
#####################################
### Key impeding matrix Seurat merging
# Function to remove DoubletFinder columns in Seurat
remove_doubletfinder_columns <- function(seurat_obj) {
  df_columns <- grep("DF.classifications|pANN", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(df_columns) > 0) {
    seurat_obj@meta.data <- seurat_obj@meta.data[, !colnames(seurat_obj@meta.data) %in% df_columns]
  }
  return(seurat_obj)
}

# Load the saved integrated batches
batch_indices <- 1:5  # Adjust based on the number of batches you have
integrated_batches <- load_integrated_batches(batch_indices)

for (i in seq_along(integrated_batches)) {
  integrated_batches[[i]] <- remove_doubletfinder_columns(integrated_batches[[i]])
}

for (i in seq_along(integrated_batches)) {
  integrated_batches[[i]] <- RenameCells(integrated_batches[[i]], add.cell.id = paste0("Batch_", i, "_"))
}


###################
#Merging
###################
merge_in_steps <- function(obj_list, step_size = 5) {
  while (length(obj_list) > 1) {
    print(paste("Current number of objects to merge:", length(obj_list)))
    merged_list <- list()
    for (i in seq(1, length(obj_list), by = step_size)) {
      batch <- obj_list[i:min(i + step_size - 1, length(obj_list))]
      print(paste("Merging batch from index", i, "to", min(i + step_size - 1, length(obj_list))))
      merged_obj <- Reduce(function(x, y) merge(x, y), batch)
      merged_list <- c(merged_list, list(merged_obj))
    }
    obj_list <- merged_list
  }
  return(obj_list[[1]])
}

combined_seurat <- merge_in_steps(integrated_batches, step_size = 5)