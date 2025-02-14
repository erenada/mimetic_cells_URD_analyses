# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
})

# Create output directories if they don't exist
dir.create("test_results/variable_genes", recursive = TRUE, showWarnings = FALSE)
dir.create("test_results/plots", recursive = TRUE, showWarnings = FALSE)

# Find the most recent test URD object
test_files <- list.files("test_data", pattern = "test_urd_object_.*\\.rds", full.names = TRUE)
if (length(test_files) == 0) {
  stop("No test URD object found. Please run create_test_urd_object.R first.")
}
latest_test_file <- test_files[which.max(file.info(test_files)$mtime)]

# Load the test URD object
message(sprintf("Loading test URD object from: %s", latest_test_file))
urd_object <- readRDS(latest_test_file)

# Calculate variable genes for all cells at once
message("\nCalculating variable genes for all cells...")
pdf("test_results/plots/variable_genes.pdf")

# Calculate variable genes using all cells
var.genes <- findVariableGenes(urd_object, 
                             cells.fit = colnames(urd_object@logupx.data), 
                             set.object.var.genes = TRUE,  # Set directly in object
                             diffCV.cutoff = 0.3,  # More permissive for test data
                             mean.min = 0.1,      
                             mean.max = 100,      
                             do.plot = TRUE)

dev.off()

message(sprintf("Variable genes found: %d", length(var.genes)))

# Save variable genes
write.table(var.genes, 
            file = "test_results/variable_genes/test_var_genes.txt",
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

# Save the updated URD object
saveRDS(urd_object, sub("\\.rds$", "_with_var_genes.rds", latest_test_file))
message("\nURD object updated with variable genes and saved") 