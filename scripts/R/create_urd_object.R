# This script is to create a URD object from the seurat object and save it as an RDS file.
# If you created the URD object before, you can skip this script.

# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
})

# Create data directory if it doesn't exist
dir.create("data", showWarnings = FALSE)

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

# Print some basic statistics about the data
message(sprintf("Number of cells: %d", ncol(seurat_object)))
message(sprintf("Number of genes: %d", nrow(seurat_object)))

# extract the counts matrix and metadata
message("Extracting counts matrix and metadata...")
counts_matrix <- GetAssayData(seurat_object, slot = "counts")
metadata <- seurat_object@meta.data

# Validate counts matrix
if (any(is.na(counts_matrix))) {
  warning("Count matrix contains NA values")
}
if (any(counts_matrix < 0)) {
  stop("Count matrix contains negative values")
}

# get the stage labels and validate
stage_labels <- as.data.frame(seurat_object@active.ident)
if (length(unique(stage_labels$`seurat_object@active.ident`)) < 2) {
  warning("Less than 2 unique stages found. URD requires multiple stages for trajectory inference.")
}

# add the stage labels to the metadata
metadata$stages <- stage_labels$`seurat_object@active.ident`

# Save intermediate files
message("Saving intermediate files...")

#save the metadata as a csv file 
write.csv(metadata, "data/metadata.csv")
message("Metadata saved.")

# save the counts matrix as a csv file
write.csv(counts_matrix, "data/counts_matrix.csv")
message("Counts matrix saved.")

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

# save the URD object with a fixed name
message("Saving URD object...")
saveRDS(urd_object, "data/initial_urd_object.rds")
message("URD object saved successfully.")

# delete the seurat object, counts matrix, and metadata from the environment to free memory
rm(seurat_object, counts_matrix, metadata)
gc() # Run garbage collection

# Verify the saved object can be loaded
message("Verifying saved URD object...")
urd_object <- readRDS("data/initial_urd_object.rds")
message("URD object verified successfully.")

# Print final summary
message("\nURD object creation complete!")
message(sprintf("Number of cells in URD object: %d", ncol(urd_object@count.data)))
message(sprintf("Number of genes in URD object: %d", nrow(urd_object@count.data)))
message(sprintf("Number of stages: %d", length(unique(urd_object@group.ids$stage))))

