# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(MASS)  # Required for fitdistr
})

# Function to safely handle PDF creation
safe_pdf <- function(filename, expr) {
  pdf(filename)
  tryCatch(
    expr,
    finally = {
      while (dev.cur() > 1) dev.off()
    }
  )
}

# Create all necessary directories
dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("results/variable_genes", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots/variable_genes", recursive = TRUE, showWarnings = FALSE)
dir.create("results/stats", recursive = TRUE, showWarnings = FALSE)

# Check if URD object exists
if (!file.exists("data/initial_urd_object.rds")) {
  stop("URD object not found. Please run create_urd_object.R first.")
}

# load the URD object
message("Loading URD object...")
urd_object <- readRDS("data/initial_urd_object.rds")

# Get stages and print summary
stages <- sort(unique(urd_object@group.ids$stage))
message("Processing stages: ", paste(stages, collapse=", "))
message(sprintf("Total number of stages: %d\n", length(stages)))

# Function to filter cells with customizable criteria
filterCells <- function(urd_obj, cells, min_counts = 100, min_genes = 50) {
  # Calculate statistics for each cell
  total_counts <- colSums(urd_obj@count.data[, cells, drop=FALSE])
  total_genes <- colSums(urd_obj@count.data[, cells, drop=FALSE] > 0)
  
  # Print before filtering statistics
  message(sprintf("Before filtering:"))
  message(sprintf("  - Number of cells: %d", length(cells)))
  message(sprintf("  - Median counts per cell: %.1f", median(total_counts)))
  message(sprintf("  - Median genes per cell: %.1f", median(total_genes)))
  message(sprintf("  - Min counts per cell: %.1f", min(total_counts)))
  message(sprintf("  - Min genes per cell: %.1f", min(total_genes)))
  
  # Check if filtering is needed
  needs_filtering <- any(total_counts < min_counts | total_genes < min_genes)
  
  if (!needs_filtering) {
    message("\nNo filtering needed - all cells already meet the minimum criteria:")
    message(sprintf("  - All cells have >= %d counts", min_counts))
    message(sprintf("  - All cells have >= %d genes", min_genes))
    
    return(list(
      cells = cells,
      stats = data.frame(
        total_cells = length(cells),
        passing_cells = length(cells),
        min_count = min(total_counts),
        median_count = median(total_counts),
        max_count = max(total_counts),
        min_genes = min(total_genes),
        median_genes = median(total_genes),
        max_genes = max(total_genes)
      )
    ))
  }
  
  # If we reach here, filtering is needed
  valid_cells <- cells[total_counts >= min_counts & total_genes >= min_genes]
  
  # Calculate after filtering statistics
  filtered_counts <- total_counts[total_counts >= min_counts & total_genes >= min_genes]
  filtered_genes <- total_genes[total_counts >= min_counts & total_genes >= min_genes]
  
  # Print after filtering statistics
  message(sprintf("\nAfter filtering (min_counts=%d, min_genes=%d):", min_counts, min_genes))
  message(sprintf("  - Number of cells: %d", length(valid_cells)))
  message(sprintf("  - Median counts per cell: %.1f", median(filtered_counts)))
  message(sprintf("  - Median genes per cell: %.1f", median(filtered_genes)))
  message(sprintf("  - Cells removed: %d (%.1f%%)", 
                 length(cells) - length(valid_cells),
                 100 * (length(cells) - length(valid_cells)) / length(cells)))
  
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
all_stage_stats <- list()  # Store statistics for all stages

# Calculate variable genes for each stage
safe_pdf("results/plots/variable_genes_by_stage.pdf", {
  for(stage in stages) {
    message(sprintf("\nProcessing stage: %s", stage))
    
    # Get and filter cells for this stage
    stage_cells <- rownames(urd_object@group.ids)[which(urd_object@group.ids$stage == stage)]
    result <- filterCells(urd_object, stage_cells, min_counts = 100, min_genes = 50)
    cells.each.stage[[stage]] <- result$cells
    
    # Store statistics with stage information
    stage_stats <- result$stats
    stage_stats$stage <- stage
    all_stage_stats[[stage]] <- stage_stats
    
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
})

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

# Save the updated URD object
saveRDS(urd_object, "data/urd_object_with_var_genes.rds")

# Combine and save all stage statistics
all_stats_df <- do.call(rbind, all_stage_stats)
write.csv(all_stats_df, 
          file = "results/stats/cell_filtering_statistics.csv",
          row.names = FALSE)

message("\nVariable genes analysis complete!")
message("Results saved to:")
message("- URD object: data/urd_object_with_var_genes.rds")
message("- Combined variable genes: results/variable_genes/all_stages_combined_var_genes.txt")
message("- Statistics: results/stats/cell_filtering_statistics.csv")
message("- Stage-specific results: results/variable_genes/")
message("- Plots: results/plots/variable_genes/")
