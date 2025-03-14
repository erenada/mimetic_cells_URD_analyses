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

# Get dataset characteristics for parameter selection
n_cells <- ncol(urd_object@logupx.data)
n_stages <- length(unique(urd_object@group.ids$stage))
cells_per_stage <- table(urd_object@group.ids$stage)
min_cells_per_stage <- min(cells_per_stage)
median_cells_per_stage <- median(cells_per_stage)
stage_size_variation <- sd(cells_per_stage) / mean(cells_per_stage)

message(sprintf("\nDataset characteristics:"))
message(sprintf("Number of cells: %d", n_cells))
message(sprintf("Number of stages: %d", n_stages))
message(sprintf("Cells per stage:"))
print(cells_per_stage)

# Calculate minimum cells flooded parameter
# Base value using square root scaling of median stage size
base_min_cells <- ceiling(sqrt(median_cells_per_stage) / 2)

# Adjust for stage size variation
variation_adjustment <- 1 + stage_size_variation
min_cells_flooded <- ceiling(base_min_cells * variation_adjustment)

# Ensure reasonable bounds
min_cells_flooded <- max(3, min(min_cells_flooded, ceiling(median_cells_per_stage * 0.05)))

message(sprintf("\nParameter selection (data-driven approach):"))
message(sprintf("Base minimum cells: %d", base_min_cells))
message(sprintf("Stage size variation adjustment: %.2f", variation_adjustment))
message(sprintf("Final minimum cells flooded: %d", min_cells_flooded))

# Identify root stage and cells
root_stage <- "Immature"  # Adjust this based on your dataset
root_cells <- which(urd_object@group.ids$stage == root_stage)

if(length(root_cells) == 0) {
    stop(sprintf("No cells found in root stage '%s'", root_stage))
}

message(sprintf("\nRoot population:"))
message(sprintf("Stage: %s", root_stage))
message(sprintf("Number of root cells: %d", length(root_cells)))

# Calculate pseudotime
message("\nCalculating pseudotime from root cells...")
urd_object <- floodPseudotime(
    object = urd_object,
    root.cells = root_cells,
    n = 100,  # Number of random walks
    minimum.cells.flooded = min_cells_flooded,
    verbose = TRUE
)

# Process flood results
message("\nProcessing flood results...")
urd_object <- floodPseudotimeProcess(
    object = urd_object,
    floods.name = "pseudotime",
    max.frac.NA = 0.4,
    pseudotime.fun = mean,
    stability.div = 20,
    verbose = TRUE
)

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
        sprintf("%d", min_cells_flooded),
        sprintf("%.2f", stage_size_variation),
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