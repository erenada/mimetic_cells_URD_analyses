# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(RColorBrewer)
})

# Create all necessary directories
dir.create("./data", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/dimensionality_reduction", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/plots/dimensionality_reduction/pca", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/plots/dimensionality_reduction/tsne", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/plots/dimensionality_reduction/outliers", recursive = TRUE, showWarnings = FALSE)

# Load the URD object with variable genes
urd_object <- readRDS("data/urd_object_with_var_genes.rds")

message("Starting dimensionality reduction analysis...")

# Get dataset characteristics
n_cells <- ncol(urd_object@logupx.data)
n_genes <- nrow(urd_object@logupx.data)
n_stages <- length(unique(urd_object@group.ids$stage))

message(sprintf("\nDataset characteristics:"))
message(sprintf("Number of cells: %d", n_cells))
message(sprintf("Number of genes: %d", n_genes))
message(sprintf("Number of stages: %d", n_stages))
message(sprintf("Number of variable genes: %d", length(urd_object@var.genes)))

# Calculate dataset complexity metrics
cells_per_stage <- table(urd_object@group.ids$stage)
min_cells_per_stage <- min(cells_per_stage)

# Adjust complexity factor calculation for larger datasets
complexity_factor <- if(n_cells > 10000) {
    (n_stages * log2(min_cells_per_stage)) / log2(n_cells)
} else {
    n_stages / log2(n_cells)
}

message("\nDataset complexity metrics:")
message(sprintf("Number of stages: %d", n_stages))
message(sprintf("Minimum cells per stage: %d", min_cells_per_stage))
message(sprintf("Complexity factor: %.2f", complexity_factor))

# 1. PCA Parameters
# -----------------------------
message("\n1. PCA Analysis")
# Try different mp.factor values
mp_factors <- c(2, 3, 4)  # Testing only 2, 3, and 4
pca_results <- list()

for(mp in mp_factors) {
  png(sprintf("results/plots/dimensionality_reduction/pca/pca_mp%.1f.png", mp),
      width = 800, height = 600, res = 100)
  urd_object <- calcPCA(urd_object, mp.factor = mp)
  pca_results[[as.character(mp)]] <- length(urd_object@pca.sig)
  pcSDPlot(urd_object)
  title(main = sprintf("mp.factor = %.1f", mp))
  dev.off()
}

# Choose optimal mp.factor
message("\nPCA significant components with different mp.factor values:")
for(mp in names(pca_results)) {
    message(sprintf("mp.factor = %s: %d significant PCs", mp, pca_results[[mp]]))
}

# Use mp.factor = 2 as default unless results suggest otherwise
chosen_mp <- 2
urd_object <- calcPCA(urd_object, mp.factor = chosen_mp)
n_sig_pcs <- length(urd_object@pca.sig)

# 2. tSNE Parameters
# -----------------------------
message("\n2. tSNE Analysis")

# Calculate perplexity based on dataset size
min_perp <- 5
max_perp <- min(50, floor(sqrt(n_cells)/3))
perplexity_values <- unique(round(c(
    min_perp,
    floor(max_perp/2),
    max_perp
)))

message(sprintf("Testing perplexity values: %s", paste(perplexity_values, collapse=", ")))

# Store results
tsne_results <- list()
successful_perp <- c()

for(perp in perplexity_values) {
    message(sprintf("\nProcessing tSNE with perplexity %d...", perp))
    
    # Set random seed for reproducibility
    set.seed(123)
    
    # Try tSNE calculation
    tsne_result <- tryCatch({
        urd_object <- calcTsne(urd_object, perplexity = perp, theta = 0.5)
        if(is.null(urd_object@tsne.y)) {
            stop("tSNE calculation produced NULL result")
        }
        TRUE
    }, error = function(e) {
        message(sprintf("Error in tSNE calculation: %s", e$message))
        return(FALSE)
    })
    
    if(tsne_result) {
        message(sprintf("tSNE calculation successful for perplexity %d", perp))
        tsne_results[[as.character(perp)]] <- urd_object@tsne.y
        successful_perp <- c(successful_perp, perp)
    } else {
        message(sprintf("tSNE calculation failed for perplexity %d", perp))
    }
    
    # Force garbage collection
    gc()
}

# Report results
message("\ntSNE Analysis Summary:")
if(length(successful_perp) > 0) {
    message(sprintf("Successful perplexity values: %s", paste(successful_perp, collapse=", ")))
} else {
    stop("No successful tSNE calculations. Cannot proceed.")
}

# Store the best tSNE result
if(length(successful_perp) > 0) {
    max_successful_perp <- max(successful_perp)
    urd_object@tsne.y <- tsne_results[[as.character(max_successful_perp)]]
    message(sprintf("\nUsing final tSNE results from perplexity %d", max_successful_perp))
} else {
    stop("No successful tSNE calculations. Cannot proceed.")
}

# Force garbage collection
gc()

# 3. Clustering Parameters
# -----------------------------
message("\n3. Graph-based Clustering")

# Calculate number of nearest neighbors based on dataset size AND complexity
base_min_nn <- ceiling(sqrt(n_cells)/2)
base_max_nn <- ceiling(sqrt(n_cells))

# Adjust nn based on complexity:
complexity_adjustment <- 1 / (1 + log2(complexity_factor))
min_nn <- ceiling(base_min_nn * complexity_adjustment)
max_nn <- ceiling(base_max_nn * complexity_adjustment)

# Ensure minimum reasonable values
min_nn <- max(min_nn, 10)  # Never go below 10 neighbors
max_nn <- max(max_nn, min_nn * 1.5)  # Ensure reasonable range

nn_values <- unique(round(c(min_nn, (min_nn + max_nn)/2, max_nn)))

message("\nNearest neighbor selection:")
message(sprintf("Base range (from sqrt rule): %d-%d", base_min_nn, base_max_nn))
message(sprintf("Complexity adjustment factor: %.2f", complexity_adjustment))
message(sprintf("Final adjusted range: %d-%d", min_nn, max_nn))
message(sprintf("Testing nearest neighbor values: %s", paste(nn_values, collapse=", ")))

# Calculate clustering for each nn value
clustering_results <- list()
for(nn in nn_values) {
  message(sprintf("\nCalculating clustering with %d nearest neighbors...", nn))
  set.seed(123)
  # Calculate clustering for this nn value
  urd_object <- graphClustering(urd_object, 
                             dim.use = "pca", 
                             num.nn = nn,
                             do.jaccard = TRUE, 
                             method = "Louvain")
  
  # Store clustering results
  cluster_name <- sprintf("Louvain-%d", nn)
  if(!is.null(urd_object@group.ids[[cluster_name]])) {
    message(sprintf("Successfully calculated clusters for nn=%d", nn))
    clustering_results[[as.character(nn)]] <- urd_object@group.ids[[cluster_name]]
  } else {
    message(sprintf("Warning: Clustering failed for nn=%d", nn))
  }
}

# Force garbage collection
gc()

# 4. kNN and Outlier Detection
# -----------------------------
message("\n4. kNN and Outlier Detection")
# Calculate optimal kNN value based on dataset size
if(n_cells < 1000) {
    nn_value <- ceiling(sqrt(n_cells))
} else if(n_cells < 10000) {
    nn_value <- ceiling(sqrt(n_cells)/2)
} else {
    nn_value <- 100
}

# Make sure we calculate enough neighbors for the 20th neighbor
nn_value <- max(nn_value, 20)
urd_object <- calcKNN(urd_object, nn = nn_value)
dist_matrix <- urd_object@knn$nn.dists

# Function to calculate data-driven bounds
calculate_data_driven_bounds <- function(dist_matrix, nn_value) {
    # Get distances - use 1st and 20th neighbors (as in URD)
    d1 <- dist_matrix[,1]   # Distance to 1st neighbor
    d20 <- dist_matrix[,20] # Distance to 20th neighbor (more distant neighborhood)
    
    # Calculate the ratio of distances and fit a line through origin
    slope <- median(d20/d1)
    
    # Calculate residuals from this line
    predicted <- slope * d1
    residuals <- d20 - predicted
    
    # Calculate the spread of residuals relative to d1
    rel_residuals <- residuals / d1
    
    # Set bounds based on the spread of relative residuals
    rel_spread <- 0.10  # Using 10% spread as it worked well
    
    # Calculate bounds that scale with d1
    slope_r <- slope * (1 + rel_spread)
    slope_b <- slope * (1 - rel_spread)
    int_r <- 0  # Intercept through origin
    int_b <- 0  # Intercept through origin
    
    # Calculate x.max using 98th percentile of d1
    x_max <- quantile(d1, 0.98)
    
    # Create diagnostic plot
    png("results/plots/dimensionality_reduction/outliers/data_driven_fit.png",
        width = 800, height = 600, res = 100)
    plot(d1, d20, 
         xlab = "Distance to neighbor 1",
         ylab = "Distance to neighbor 20",
         main = "Data-Driven Boundary Fitting\n(URD-style: 1st vs 20th neighbor)")
    abline(0, slope, col = "black", lty = 2)  # Base fit
    abline(int_r, slope_r, col = "red")  # Upper bound
    abline(int_b, slope_b, col = "blue")  # Lower bound
    abline(v = x_max, col = "green")  # x.max
    dev.off()
    
    # Print fit statistics
    message("\nData-driven bound statistics:")
    message(sprintf("Median d20/d1 ratio (base slope): %.3f", slope))
    message(sprintf("Upper bound slope: %.3f", slope_r))
    message(sprintf("Lower bound slope: %.3f", slope_b))
    message(sprintf("x.max (98th percentile): %.2f", x_max))
    
    return(list(
        slope.r = slope_r,
        int.r = int_r,
        slope.b = slope_b,
        int.b = int_b,
        x.max = x_max,
        base_slope = slope
    ))
}

# Function to identify outliers
identify_outliers <- function(d1, d20, params) {
    # A cell is an outlier if:
    # 1. d1 > x.max OR
    # 2. d20 > slope.r * d1 + int.r OR
    # 3. d20 < slope.b * d1 + int.b
    
    outliers <- which(
        d1 > params$x.max |  # First neighbor too far
        d20 > (params$slope.r * d1 + params$int.r) |  # Above upper bound
        d20 < (params$slope.b * d1 + params$int.b)    # Below lower bound
    )
    
    return(outliers)
} 

# Calculate data-driven parameters
bounds <- calculate_data_driven_bounds(dist_matrix, nn_value)

# Create parameter combinations around the data-driven values
parameter_combinations <- list(
    # Single parameter set with 10% spread and 98th percentile x.max
    list(slope.r = bounds$base_slope * 1.10, 
         int.r = bounds$int.r, 
         slope.b = bounds$base_slope * 0.90, 
         int.b = bounds$int.b,
         x.max = bounds$x.max,
         description = "URD-style bounds (10% spread, 98th percentile x.max)")
)

outliers_results <- list()
for(i in seq_along(parameter_combinations)) {
    params <- parameter_combinations[[i]]
    
    # Get the distances
    d1 <- dist_matrix[,1]
    d20 <- dist_matrix[,20]
    
    # Identify outliers using our custom function
    outlier_indices <- identify_outliers(d1, d20, params)
    outliers <- colnames(urd_object@logupx.data)[outlier_indices]
    
    # Create the plot
    png(sprintf("results/plots/dimensionality_reduction/outliers/outliers_set%d.png", i),
        width = 800, height = 600, res = 100)
    plot(d1, d20,
         xlab = "Distance to neighbor 1",
         ylab = "Distance to neighbor 20",
         main = sprintf("%s\n(Parameter Set %d)", params$description, i))
    points(d1[outlier_indices], d20[outlier_indices], col = "red", pch = 16)
    abline(params$int.r, params$slope.r, col = "red")  # Upper bound
    abline(params$int.b, params$slope.b, col = "blue")  # Lower bound
    abline(v = params$x.max, col = "green")  # x.max
    dev.off()
    
    outliers_results[[i]] <- list(
        params = params,
        outliers = outliers,
        percent = 100 * length(outliers) / n_cells,
        description = params$description
    )
    
    # Print detailed results for each parameter set
    message(sprintf("\nParameter Set %d (%s):", i, params$description))
    message(sprintf("Slope (red): %.3f, Intercept (red): %.3f", params$slope.r, params$int.r))
    message(sprintf("Slope (blue): %.3f, Intercept (blue): %.3f", params$slope.b, params$int.b))
    message(sprintf("Outliers: %d (%.2f%%)", length(outliers), 100 * length(outliers) / n_cells))
}

# Modify the best result selection criteria for larger datasets
best_result <- NULL
target_outlier_range <- if(n_cells > 10000) {
    c(1, 20)  # Allow up to 20% outliers for larger datasets
} else {
    c(1, 10)  # Original range for smaller datasets
}

# Compare results and select best parameters
for(i in seq_along(outliers_results)) {
    if(outliers_results[[i]]$percent >= target_outlier_range[1] && 
       outliers_results[[i]]$percent <= target_outlier_range[2]) {
        best_result <- outliers_results[[i]]
        message(sprintf("\nSelected Parameter Set %d (%s) as it gives reasonable outlier percentage.", 
                      i, outliers_results[[i]]$description))
        break
    }
}

# If no good result found, use most lenient set
if(is.null(best_result)) {
    best_result <- outliers_results[[length(outliers_results)]]
    message(sprintf("\nNo parameter set gave %.1f-%.1f%% outliers. Using most lenient set (%s).", 
                   target_outlier_range[1], target_outlier_range[2],
                   best_result$description))
}

best_outliers <- best_result$outliers

# Analyze outlier distribution across stages
message("\nAnalyzing outlier distribution across stages:")
stage_counts <- table(urd_object@group.ids$stage)
outlier_stages <- urd_object@group.ids$stage[colnames(urd_object@logupx.data) %in% best_outliers]
outlier_by_stage <- table(outlier_stages)
percent_by_stage <- (outlier_by_stage / stage_counts[names(outlier_by_stage)]) * 100

# Create a summary data frame
stage_summary <- data.frame(
    Stage = names(stage_counts),
    Total_Cells = as.numeric(stage_counts),
    Outlier_Cells = as.numeric(outlier_by_stage[names(stage_counts)]),
    Outlier_Percentage = as.numeric(percent_by_stage[names(stage_counts)])
)
stage_summary$Outlier_Cells[is.na(stage_summary$Outlier_Cells)] <- 0
stage_summary$Outlier_Percentage[is.na(stage_summary$Outlier_Percentage)] <- 0

# Print stage-wise statistics
message("\nOutlier distribution by stage:")
for(i in 1:nrow(stage_summary)) {
    message(sprintf("%s: %d outliers out of %d cells (%.2f%%)", 
                   stage_summary$Stage[i],
                   stage_summary$Outlier_Cells[i],
                   stage_summary$Total_Cells[i],
                   stage_summary$Outlier_Percentage[i]))
}

# Save stage statistics to file
write.csv(stage_summary, 
          "results/dimensionality_reduction/outlier_stage_statistics.csv", 
          row.names = FALSE, 
          quote = FALSE)

# Save outlier information
writeLines(c(
    sprintf("# Outlier Detection Parameters Used:"),
    sprintf("# %s", best_result$description),
    sprintf("# Slope (red): %.3f, Intercept (red): %.3f", best_result$params$slope.r, best_result$params$int.r),
    sprintf("# Slope (blue): %.3f, Intercept (blue): %.3f", best_result$params$slope.b, best_result$params$int.b),
    sprintf("# Outlier percentage: %.2f%%", best_result$percent),
    sprintf("# Total cells analyzed: %d", n_cells),
    sprintf("# Number of outlier cells: %d", length(best_outliers)),
    "",
    "# Stage-wise Outlier Statistics:",
    sprintf("# %s: %d outliers (%.2f%%)", 
            stage_summary$Stage,
            stage_summary$Outlier_Cells,
            stage_summary$Outlier_Percentage),
    "",
    "# Dataset Statistics:",
    sprintf("# Minimum distance: %.2f", min(dist_matrix)),
    sprintf("# Median distance: %.2f", median(dist_matrix)),
    sprintf("# 95th percentile (x.max): %.2f", bounds$x.max),
    sprintf("# Maximum distance: %.2f", max(dist_matrix)),
    "",
    "# Cell IDs identified as outliers:",
    best_outliers
), "results/dimensionality_reduction/outlier_detection_parameters.txt")

# Also save just the outlier cell IDs in a separate file
write.table(best_outliers, 
           file = "results/dimensionality_reduction/outlier_cells.txt", 
           row.names = FALSE, 
           col.names = FALSE, 
           quote = FALSE)

message(sprintf("\nIdentified %d outliers (%.1f%%) using %s", 
                length(best_outliers), 
                100 * length(best_outliers) / n_cells,
                best_result$description))

# Create clean subset
cells.keep <- setdiff(colnames(urd_object@logupx.data), best_outliers)
urd_object_clean <- urdSubset(urd_object, cells.keep = cells.keep)

# Restore clustering results
for(nn in names(clustering_results)) {
    cluster_name <- sprintf("Louvain-%d", as.numeric(nn))
    urd_object@group.ids[[cluster_name]] <- clustering_results[[nn]]
    urd_object_clean@group.ids[[cluster_name]] <- clustering_results[[nn]][cells.keep]
}

# Save results
saveRDS(urd_object, "data/urd_object_with_dimred.rds")
saveRDS(urd_object_clean, "data/urd_object_clean.rds")

# Save parameter choices
parameter_summary <- data.frame(
    parameter = c("PCA mp.factor", 
                 "Significant PCs",
                 "tSNE perplexity range",
                 "Clustering nn range",
                 "kNN value",
                 "Outlier percentage"),
    value = c(sprintf("%.1f", chosen_mp),
             sprintf("%d", n_sig_pcs),
             sprintf("%d-%d", min(perplexity_values), max(perplexity_values)),
             sprintf("%d-%d", min(nn_values), max(nn_values)),
             sprintf("%d", nn_value),
             sprintf("%.1f%%", 100 * length(best_outliers) / n_cells))
)
write.csv(parameter_summary, 
          "results/dimensionality_reduction/parameters.csv", 
          row.names = FALSE, 
          quote = FALSE)

message("\nDimensionality reduction analysis complete!")
message("Results saved to:")
message("- URD object with dimensionality reduction: data/urd_object_with_dimred.rds")
message("- Clean URD object: data/urd_object_clean.rds")
message("- Parameter summary: results/dimensionality_reduction/parameters.csv")
message("- Plots: results/plots/dimensionality_reduction/") 

