# This script creates a small test URD object by subsampling the original data
# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
})

# Create test data directory if it doesn't exist
dir.create("test_data", showWarnings = FALSE)
dir.create("test_results", showWarnings = FALSE)

# Function to check if file exists
check_file_exists <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Input file not found: %s", filepath))
  }
}

# Input file path
input_file <- "data/meclo_ht2-ht6_seurat_table_20240503_labelled-ExportForEren20240718.rds"
check_file_exists(input_file)

# load seurat data with error handling
tryCatch({
  message("Loading Seurat object...")
  seurat_object <- readRDS(input_file)
  
  # Validate that it's a Seurat object
  if (!inherits(seurat_object, "Seurat")) {
    stop("Input object is not a Seurat object")
  }
}, error = function(e) {
  stop(sprintf("Error loading Seurat object: %s", e$message))
})

# Print original data statistics
message("Original dataset characteristics:")
message(sprintf("Number of cells: %d", ncol(seurat_object)))
message(sprintf("Number of genes: %d", nrow(seurat_object)))

# Subsample cells (10% of original data or 1000 cells, whichever is smaller)
set.seed(42)  # for reproducibility
n_cells_to_sample <- min(1000, round(ncol(seurat_object) * 0.1))
cells_to_keep <- sample(colnames(seurat_object), n_cells_to_sample)
seurat_object_small <- subset(seurat_object, cells = cells_to_keep)

# Print subsampled data statistics
message("\nSubsampled dataset characteristics:")
message(sprintf("Number of cells: %d", ncol(seurat_object_small)))
message(sprintf("Number of genes: %d", nrow(seurat_object_small)))

# extract the counts matrix and metadata
message("Extracting counts matrix and metadata...")
counts_matrix <- GetAssayData(seurat_object_small, slot = "counts")
metadata <- seurat_object_small@meta.data

# Validate counts matrix
if (any(is.na(counts_matrix))) {
  warning("Count matrix contains NA values")
}
if (any(counts_matrix < 0)) {
  stop("Count matrix contains negative values")
}

# get the stage labels and validate
stage_labels <- as.data.frame(seurat_object_small@active.ident)
if (length(unique(stage_labels$`seurat_object_small@active.ident`)) < 2) {
  warning("Less than 2 unique stages found. URD requires multiple stages for trajectory inference.")
}

# add the stage labels to the metadata
metadata$stages <- stage_labels$`seurat_object_small@active.ident`

# Save intermediate files with timestamps
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

# create a URD object with error handling
message("Creating URD object...")
tryCatch({
  urd_object <- createURD(counts_matrix, metadata)
}, error = function(e) {
  stop(sprintf("Error creating URD object: %s", e$message))
})

## Copy stage from @meta to @group.ids
message("Processing stage labels...")
urd_object@group.ids$stage <- as.character(urd_object@meta[rownames(urd_object@group.ids),"stages"])

# remove the spaces in the stage labels with underscores
urd_object@group.ids$stage <- gsub(" ", "_", urd_object@group.ids$stage)

# Print summary of stages
message("Stage summary:")
print(table(urd_object@group.ids$stage))

# save the URD object
message("Saving URD object...")
saveRDS(urd_object, sprintf("test_data/test_urd_object_%s.rds", timestamp))
message("URD object saved successfully.")

# delete objects to free memory
rm(seurat_object, seurat_object_small, counts_matrix, metadata)
gc() # Run garbage collection

# Verify the saved object can be loaded
message("Verifying saved URD object...")
test_urd_object <- readRDS(sprintf("test_data/test_urd_object_%s.rds", timestamp))
message("URD object verified successfully.")

# Print final summary
message("\nTest URD object creation complete!")
message(sprintf("Number of cells in URD object: %d", ncol(test_urd_object@count.data)))
message(sprintf("Number of genes in URD object: %d", nrow(test_urd_object@count.data)))
message(sprintf("Number of stages: %d", length(unique(test_urd_object@group.ids$stage)))) 