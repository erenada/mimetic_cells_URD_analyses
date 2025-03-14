# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(RColorBrewer)
})

# Create all necessary directories
dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("results/pseudotime", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots/pseudotime", recursive = TRUE, showWarnings = FALSE)

# Load the URD object with diffusion map
urd_object <- readRDS("data/urd_object_with_dm.rds")

message("Starting pseudotime calculation...")

# Verify diffusion map
if (is.null(urd_object@dm)) {
    stop("No diffusion map found in URD object. Please run diffusion map calculation first.")
}

if (is.null(urd_object@dm@eigenvectors) || ncol(urd_object@dm@eigenvectors) < 2) {
    stop("Invalid diffusion map: eigenvectors not properly calculated")
}

# Get dataset characteristics for parameter selection
n_cells <- ncol(urd_object@logupx.data)
n_stages <- length(unique(urd_object@group.ids$stage))
cells_per_stage <- table(urd_object@group.ids$stage)

message(sprintf("\nDataset characteristics:"))
message(sprintf("Number of cells: %d", n_cells))
message(sprintf("Number of stages: %d", n_stages))
message(sprintf("Cells per stage:"))
print(cells_per_stage)

# Identify root stage and cells
root_stage <- "Immature"  # Adjust this based on your dataset

# Check if the stage exists
if (!(root_stage %in% urd_object@group.ids$stage)) {
    stop(sprintf("Root stage '%s' not found in data. Available stages: %s", 
                 root_stage, 
                 paste(unique(urd_object@group.ids$stage), collapse=", ")))
}

# Get cell names instead of indices for root cells
root_cells <- rownames(urd_object@group.ids)[urd_object@group.ids$stage == root_stage]

if(length(root_cells) == 0) {
    stop(sprintf("No cells found in root stage '%s'", root_stage))
}

message(sprintf("\nRoot population:"))
message(sprintf("Stage: %s", root_stage))
message(sprintf("Number of root cells: %d", length(root_cells)))

# Calculate pseudotime with error handling
message("\nCalculating pseudotime from root cells...")
tryCatch({
    # Run flood simulations and save results
    flood_results <- floodPseudotime(
        object = urd_object,
        root.cells = root_cells,  # Now using cell names instead of indices
        n = 100,  # Number of random walks as per URD documentation
        minimum.cells.flooded = 2,  # Set to 2 as per URD documentation
        verbose = TRUE
    )
    
    # Save the flood results
    saveRDS(flood_results, "results/pseudotime/flood_results.rds")
    message("Flood results saved to: results/pseudotime/flood_results.rds")
    
}, error = function(e) {
    stop(sprintf("Error in floodPseudotime: %s", e$message))
})

# Process flood results
message("\nProcessing flood results...")
tryCatch({
    # Process the floods to generate pseudotime
    urd_object <- floodPseudotimeProcess(
        object = urd_object,
        floods = flood_results,  # Use the saved flood results
        floods.name = "pseudotime",
        max.frac.NA = 0.4,     # Default value from documentation
        pseudotime.fun = mean,  # Default and only validated function
        stability.div = 10,     # Default value from documentation
        verbose = TRUE
    )
}, error = function(e) {
    stop(sprintf("Error in floodPseudotimeProcess: %s", e$message))
})

# Check pseudotime calculation success
if (is.null(urd_object@pseudotime)) {
    stop("Pseudotime calculation failed - no pseudotime values found in URD object")
}

# Calculate stage statistics
stage_stats <- tapply(urd_object@pseudotime$pseudotime, urd_object@group.ids$stage, 
                     function(x) c(mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE)))
stage_stats_df <- do.call(rbind, stage_stats)
stage_order <- order(stage_stats_df[,1])

message("\nPseudotime statistics by stage:")
print(stage_stats_df[stage_order,])

# Create visualizations
message("\nGenerating visualization plots...")

# 1. Stage-wise pseudotime distribution
png("results/plots/pseudotime/stage_pseudotime_distribution.png",
    width = 1200, height = 800, res = 150)
par(mar = c(10, 4, 4, 2))
boxplot(pseudotime ~ stage,
        data = data.frame(
            pseudotime = urd_object@pseudotime$pseudotime,
            stage = urd_object@group.ids$stage
        ),
        las = 2,
        main = "Pseudotime Distribution by Stage",
        ylab = "Pseudotime",
        col = brewer.pal(8, "Set2"))
dev.off()

# 2. Pseudotime on diffusion components
png("results/plots/pseudotime/pseudotime_diffusion_components.png",
    width = 1200, height = 1200, res = 150)
plotDimArray(
    object = urd_object,
    reduction.use = "dm",
    dims.to.plot = 1:16,
    label = "pseudotime",
    plot.title = "",
    outer.title = "Pseudotime on Diffusion Components",
    legend = TRUE,
    alpha = 0.6
)
dev.off()

# 3. Stability assessment
png("results/plots/pseudotime/pseudotime_stability.png",
    width = 800, height = 800, res = 150)
plotStabilityOverall(urd_object)
dev.off()

# Save results
message("\nSaving results...")

# Save the URD object with pseudotime
saveRDS(urd_object, "data/urd_object_with_pseudotime.rds")

# Save parameter summary
parameter_summary <- data.frame(
    parameter = c(
        "Dataset size",
        "Number of stages",
        "Root stage",
        "Root population size",
        "Minimum cells flooded",
        "Stage size variation",
        "Number of simulations"
    ),
    value = c(
        sprintf("%d cells", n_cells),
        sprintf("%d", n_stages),
        root_stage,
        sprintf("%d cells", length(root_cells)),
        "2",
        "0.00",
        "100"
    )
)
write.csv(parameter_summary,
          "results/pseudotime/parameters.csv",
          row.names = FALSE,
          quote = FALSE)

# Save stage statistics
write.csv(stage_stats_df,
          "results/pseudotime/stage_statistics.csv",
          row.names = TRUE)

message("\nPseudotime analysis complete!")
message("Results saved to:")
message("- URD object: data/urd_object_with_pseudotime.rds")
message("- Parameter summary: results/pseudotime/parameters.csv")
message("- Stage statistics: results/pseudotime/stage_statistics.csv")
message("- Plots: results/plots/pseudotime/") 