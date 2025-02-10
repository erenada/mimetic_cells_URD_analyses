# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(MASS)  # Required for fitdistr
})

# load the URD object
urd_object <- readRDS("data/initial_urd_object_20250210_1206.rds")

# Create output directories if they don't exist
dir.create("results/variable_genes", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

# Get stages and print summary
stages <- sort(unique(urd_object@group.ids$stage))
message("Processing stages: ", paste(stages, collapse=", "))
message(sprintf("Total number of stages: %d\n", length(stages)))

# Function to filter cells with customizable criteria
filterCells <- function(urd_obj, cells, min_counts = 100, min_genes = 50) {
  # Calculate statistics for each cell
  total_counts <- colSums(urd_obj@count.data[, cells, drop=FALSE])
  total_genes <- colSums(urd_obj@count.data[, cells, drop=FALSE] > 0)
  
  # Apply filters
  valid_cells <- cells[total_counts >= min_counts & total_genes >= min_genes]
  
  # Return both valid cells and filtering statistics
  return(list(
    cells = valid_cells,
    stats = data.frame(
      total_cells = length(cells),
      passing_cells = length(valid_cells),
      min_count = min(total_counts),
      median_count = median(total_counts),
      max_count = max(total_counts),
      min_genes = min(total_genes),
      median_genes = median(total_genes),
      max_genes = max(total_genes)
    )
  ))
}

# Process each stage
message("Processing stages and calculating variable genes...")
cells.each.stage <- list()
var.genes.by.stage <- list()

# Set up PDF for all plots
pdf("results/plots/variable_genes_by_stage.pdf")

for(stage in stages) {
  message(sprintf("\nProcessing stage: %s", stage))
  
  # Get and filter cells for this stage
  stage_cells <- rownames(urd_object@group.ids)[which(urd_object@group.ids$stage == stage)]
  result <- filterCells(urd_object, stage_cells, min_counts = 100, min_genes = 50)
  cells.each.stage[[stage]] <- result$cells
  
  message(sprintf("Cells: %d (after filtering)", length(result$cells)))
  
  # Calculate variable genes
  var.genes <- findVariableGenes(urd_object, 
                               cells.fit = result$cells, 
                               set.object.var.genes = FALSE, 
                               diffCV.cutoff = 0.3,  
                               mean.min = 0.1,      
                               mean.max = 100,      
                               main.use = stage,
                               do.plot = TRUE)
  
  var.genes.by.stage[[stage]] <- var.genes
  message(sprintf("Variable genes found: %d", length(var.genes)))
  
  # Save variable genes for this stage
  write.table(var.genes, 
              file = sprintf("results/variable_genes/%s_var_genes.txt", stage),
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
}

# Close the PDF device
dev.off()

# Combine variable genes from all stages
var.genes <- sort(unique(unlist(var.genes.by.stage)))
urd_object@var.genes <- var.genes

# Save combined variable genes
write.table(var.genes, 
            file = "results/variable_genes/all_stages_combined_var_genes.txt",
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

# Print summary
message("\nAnalysis Summary:")
message(sprintf("Total variable genes: %d", length(var.genes)))
message("\nVariable genes per stage:")
for(stage in names(var.genes.by.stage)) {
  message(sprintf("%s: %d", stage, length(var.genes.by.stage[[stage]])))
}

# Save the variable genes statistics
var_genes_stats <- data.frame(
  stage = names(var.genes.by.stage),
  num_var_genes = sapply(var.genes.by.stage, length)
)
write.csv(var_genes_stats, 
          file = "results/variable_genes/variable_genes_statistics.csv",
          row.names = FALSE)

# Overwrite the original URD object
saveRDS(urd_object, "data/initial_urd_object_20250210_1206.rds")
message("\nURD object updated with variable genes")
