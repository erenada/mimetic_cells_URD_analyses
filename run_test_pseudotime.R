# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(RColorBrewer)
})

# Find the most recent test URD object with diffusion map
test_files <- list.files("test_data", pattern = "test_urd_object_with_dm\\.rds$", full.names = TRUE)
if (length(test_files) == 0) {
  stop("No URD object with diffusion map found. Please run run_test_diffusion_map.R first.")
}
latest_test_file <- test_files[which.max(file.info(test_files)$mtime)]

# Load the URD object
message(sprintf("Loading URD object from: %s", latest_test_file))
urd_object <- readRDS(latest_test_file)

# Create output directories
dir.create("test_results/plots/pseudotime", recursive = TRUE, showWarnings = FALSE)

message("Starting pseudotime calculation...")

# Get dataset characteristics
n_cells <- ncol(urd_object@logupx.data)
n_stages <- length(unique(urd_object@group.ids$stage))
cells_per_stage <- table(urd_object@group.ids$stage)

message(sprintf("\nDataset characteristics:"))
message(sprintf("Number of cells: %d", n_cells))
message(sprintf("Number of stages: %d", n_stages))
message(sprintf("Cells per stage:"))
print(cells_per_stage)

# Print all stage names for debugging
message("\nAvailable stages:")
print(levels(factor(urd_object@group.ids$stage)))

# Identify root cells (Immature stage) - direct matching
root_stage <- "Immature"  # Exact match since we can see it in the output
root_cells <- colnames(urd_object@logupx.data)[urd_object@group.ids$stage == root_stage]

if (length(root_cells) == 0) {
  stop("Could not find cells from '", root_stage, "' stage. Available stages are: ", 
       paste(sort(unique(urd_object@group.ids$stage)), collapse=", "))
}

message(sprintf("\nFound %d root cells from stage: '%s'", length(root_cells), root_stage))

# Function to determine optimal parameters based on dataset characteristics
determine_optimal_params <- function(n_cells, n_stages) {
  # Calculate dataset complexity metrics
  min_cells_per_stage <- min(cells_per_stage)
  median_cells_per_stage <- median(as.numeric(cells_per_stage))
  stage_size_variation <- sd(as.numeric(cells_per_stage)) / median_cells_per_stage
  
  # Calculate minimum cells flooded (minimum.cells.flooded parameter)
  # We want enough cells to get stable measurements but not so many that we:
  # 1. Miss small but important cell populations
  # 2. Over-smooth the trajectory
  # 3. Create computational bottlenecks
  
  # Start with square root scaling for stability
  base_min_cells <- ceiling(sqrt(median_cells_per_stage))
  
  # Adjust based on stage size variation
  # More variation = need more cells for stable measurement
  variation_adjustment <- 1 + stage_size_variation
  min_cells_flooded <- ceiling(base_min_cells * variation_adjustment)
  
  # Ensure reasonable bounds
  # Lower bound: Need at least 3 cells for stable statistics
  # Upper bound: Don't want to use more than 5% of typical stage
  min_cells_flooded <- max(3, min(
    min_cells_flooded,
    ceiling(median_cells_per_stage * 0.05)
  ))
  
  # Print detailed parameter selection rationale
  message("\nParameter selection rationale:")
  message(sprintf("Dataset complexity metrics:"))
  message(sprintf("- Minimum cells per stage: %d", min_cells_per_stage))
  message(sprintf("- Median cells per stage: %.1f", median_cells_per_stage))
  message(sprintf("- Stage size variation: %.2f", stage_size_variation))
  
  message(sprintf("\nMinimum cells flooded calculation:"))
  message(sprintf("- Base value (sqrt of median stage): %.1f", base_min_cells))
  message(sprintf("- Variation adjustment: %.2fx", variation_adjustment))
  message(sprintf("- Upper bound (5%% of median stage): %.1f", median_cells_per_stage * 0.05))
  message(sprintf("- Final minimum cells flooded: %d", min_cells_flooded))
  
  return(list(
    cells_per_waypoint = min_cells_flooded,  # Keep parameter name for compatibility
    complexity_metrics = list(
      min_cells_per_stage = min_cells_per_stage,
      median_cells_per_stage = median_cells_per_stage,
      stage_size_variation = stage_size_variation
    )
  ))
}

# Get optimal parameters
params <- determine_optimal_params(n_cells, n_stages)
message(sprintf("\nCalculated optimal parameters:"))
message(sprintf("Minimum cells flooded: %d", params$cells_per_waypoint))

# Calculate pseudotime
message("\nCalculating pseudotime from root cells...")
set.seed(123) # For reproducibility

# Flood simulations with correct parameters
message("Running flood simulations...")
flood_results <- floodPseudotime(
  object = urd_object,
  root.cells = root_cells,
  n = 100,  # Number of simulations
  minimum.cells.flooded = params$cells_per_waypoint,  # Minimum cells per waypoint
  verbose = TRUE
)

# Process the floods
message("\nProcessing flood results...")
urd_object <- floodPseudotimeProcess(
  urd_object, 
  flood_results, 
  floods.name = "pseudotime", 
  max.frac.NA = 0.4,
  pseudotime.fun = mean, 
  stability.div = 20
)

# Check pseudotime stability
message("\nChecking pseudotime stability...")
png("test_results/plots/pseudotime/stability_plot.png",
    width = 800, height = 600, res = 100)
pseudotimePlotStabilityOverall(urd_object)
dev.off()

# Save intermediate result
saveRDS(urd_object, "test_data/test_urd_object_with_pseudotime.rds")

# Create visualization plots
message("\nGenerating visualization plots...")

# 1. Pseudotime distribution plot
png("test_results/plots/pseudotime/pseudotime_distribution.png",
    width = 800, height = 600, res = 100)
plotDists(urd_object, "pseudotime", "stage", plot.title="Pseudotime by stage")
dev.off()

# 2. Pseudotime on diffusion components
png("test_results/plots/pseudotime/pseudotime_dm_components.png",
    width = 1200, height = 1200, res = 150)
plotDimArray(urd_object, 
            reduction.use = "dm", 
            dims.to.plot = 1:9,
            label = "pseudotime", 
            plot.title = "", 
            outer.title = "Pseudotime on Diffusion Components",
            legend = TRUE)
dev.off()

# Save parameter summary
parameter_summary <- data.frame(
    parameter = c("Dataset size", 
                 "Number of stages",
                 "Root cells",
                 "Minimum cells flooded",
                 "Number of simulations",
                 "Mean pseudotime",
                 "Median pseudotime",
                 "Pseudotime range"),
    value = c(sprintf("%d cells", n_cells),
              sprintf("%d", n_stages),
              sprintf("%d (Immature stage)", length(root_cells)),
              sprintf("%d", params$cells_per_waypoint),
              "100",
              sprintf("%.2f", mean(urd_object@pseudotime$pseudotime)),
              sprintf("%.2f", median(urd_object@pseudotime$pseudotime)),
              sprintf("%.2f - %.2f", min(urd_object@pseudotime$pseudotime), max(urd_object@pseudotime$pseudotime)))
)
write.csv(parameter_summary, 
          "test_results/pseudotime_parameters.csv", 
          row.names = FALSE, 
          quote = FALSE)

message("\nPseudotime calculation complete!")
message("Results saved to 'test_data/test_urd_object_with_pseudotime.rds'")
message("Parameter summary saved to 'test_results/pseudotime_parameters.csv'")
message("Plots saved in 'test_results/plots/pseudotime/'")

# Print summary statistics
message("\nPseudotime Summary Statistics:")
pseudotime_stats <- summary(urd_object@pseudotime$pseudotime)
print(pseudotime_stats)

# Print stage-wise pseudotime averages
stage_pt_means <- tapply(urd_object@pseudotime$pseudotime, urd_object@group.ids$stage, mean, na.rm=TRUE)
message("\nMean pseudotime by stage (sorted):")
print(sort(stage_pt_means)) 