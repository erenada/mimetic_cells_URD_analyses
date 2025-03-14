# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(RColorBrewer)
})

# Find the most recent test URD object
test_files <- list.files("./test/test_data", pattern = "test_urd_object_.*_with_var_genes\\.rds$", full.names = TRUE)
if (length(test_files) == 0) {
  stop("No test URD object with variable genes found. Please run run_test_linage_analysis.R first.")
}
latest_test_file <- test_files[which.max(file.info(test_files)$mtime)]

# Load the test URD object
message(sprintf("Loading test URD object from: %s", latest_test_file))
urd_object <- readRDS(latest_test_file)


# Create all necessary directories
dir.create("./test/test_data", recursive = TRUE, showWarnings = FALSE)
dir.create("./test/test_results/dimensionality_reduction", recursive = TRUE, showWarnings = FALSE)
dir.create("./test/test_results/plots/dimensionality_reduction/pca", recursive = TRUE, showWarnings = FALSE)
dir.create("./test/test_results/plots/dimensionality_reduction/tsne", recursive = TRUE, showWarnings = FALSE)
dir.create("./test/test_results/plots/dimensionality_reduction/outliers", recursive = TRUE, showWarnings = FALSE)

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

# 1. PCA Parameters
# -----------------------------
message("\n1. PCA Analysis")
# Try different mp.factor values (adjusted for test dataset)
mp_factors <- c(2, 3, 4)  # Matching main script
pca_results <- list()

for(mp in mp_factors) {
  png(sprintf("./test/test_results/plots/dimensionality_reduction/pca/pca_mp%.1f.png", mp),
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

# Calculate perplexity based on dataset size (adjusted for smaller dataset)
min_perp <- 5
max_perp <- min(30, floor(sqrt(n_cells)/3))
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

# Calculate dataset complexity metrics (adjusted for smaller dataset)
n_stages <- length(unique(urd_object@group.ids$stage))
cells_per_stage <- table(urd_object@group.ids$stage)
min_cells_per_stage <- min(cells_per_stage)
complexity_factor <- n_stages / log2(n_cells)

message("\nDataset complexity metrics:")
message(sprintf("Number of stages: %d", n_stages))
message(sprintf("Minimum cells per stage: %d", min_cells_per_stage))
message(sprintf("Complexity factor: %.2f", complexity_factor))

# Calculate number of nearest neighbors (adjusted for smaller dataset)
base_min_nn <- ceiling(sqrt(n_cells)/3)  # Reduced from /2 for smaller dataset
base_max_nn <- ceiling(sqrt(n_cells)/2)  # Reduced from sqrt(n_cells) for smaller dataset

# Adjust nn based on complexity
complexity_adjustment <- 1 / (1 + log2(complexity_factor))
min_nn <- as.integer(ceiling(base_min_nn * complexity_adjustment))
max_nn <- as.integer(ceiling(base_max_nn * complexity_adjustment))

# Ensure minimum reasonable values
min_nn <- as.integer(max(min_nn, 5))  # Reduced from 10 for smaller dataset
max_nn <- as.integer(max(max_nn, min_nn * 1.5))

# Convert to integer for unique and round operations
nn_values <- as.integer(unique(round(c(min_nn, (min_nn + max_nn)/2, max_nn))))

message("\nNearest neighbor selection:")
message(sprintf("Base range (from sqrt rule): %d-%d", as.integer(base_min_nn), as.integer(base_max_nn)))
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
# Calculate optimal kNN value (adjusted for smaller dataset)
if(n_cells < 1000) {
    nn_value <- ceiling(sqrt(n_cells))
} else if(n_cells < 10000) {
    nn_value <- ceiling(sqrt(n_cells)/2)
} else {
    nn_value <- 100
}

# Calculate kNN only once
urd_object <- calcKNN(urd_object, nn = nn_value)
dist_matrix <- urd_object@knn$nn.dists
all_distances <- as.vector(dist_matrix)
x_max_value <- quantile(all_distances, 0.95)

message(sprintf("\nDistance distribution statistics:"))
message(sprintf("Min distance: %.2f", min(all_distances)))
message(sprintf("Median distance: %.2f", median(all_distances)))
message(sprintf("95th percentile (x.max): %.2f", x_max_value))
message(sprintf("Max distance: %.2f", max(all_distances)))

message(sprintf("\nChosen parameters for kNN:"))
message(sprintf("nn_value: %d (based on dataset size)", nn_value))
message(sprintf("nn.2: %d (1/3 of nn_value)", round(nn_value/3)))
message(sprintf("x.max: %.2f (95th percentile of distances)", x_max_value))

# Create directory for outlier plots
dir.create("./test/test_results/plots/dimensionality_reduction/outliers", recursive = TRUE, showWarnings = FALSE)

# Function to calculate data-driven bounds
calculate_data_driven_bounds <- function(dist_matrix, nn_value) {
    # Get distances
    d1 <- dist_matrix[,1]  # Distance to 1st neighbor
    dn <- dist_matrix[,round(nn_value/3)]  # Distance to n/3 neighbor
    
    # Fit linear model
    lm_fit <- lm(dn ~ d1)
    
    # Get base slope and intercept
    base_slope <- coef(lm_fit)[2]
    base_intercept <- coef(lm_fit)[1]
    
    # Calculate residuals and their standard deviation
    residuals <- residuals(lm_fit)
    resid_sd <- sd(residuals)
    
    # Calculate bounds using residual distribution
    # Upper bound: 2.5 standard deviations above regression line
    slope_r <- base_slope
    int_r <- base_intercept + 2.5 * resid_sd
    
    # Lower bound: 2.5 standard deviations below regression line
    slope_b <- base_slope
    int_b <- base_intercept - 2.5 * resid_sd
    
    # Calculate x.max using robust statistics
    x_max <- quantile(d1, 0.95)
    
    # Create diagnostic plot
    plot_file <- file.path("test", "test_results", "plots", "dimensionality_reduction", "outliers", "data_driven_fit.png")
    png(plot_file, width = 800, height = 600, res = 100)
    plot(d1, dn, 
         xlab = "Distance to neighbor 1",
         ylab = sprintf("Distance to neighbor %d", round(nn_value/3)),
         main = "Data-Driven Boundary Fitting")
    abline(lm_fit, col = "black", lty = 2)  # Base fit
    abline(a = int_r, b = slope_r, col = "red")  # Upper bound
    abline(a = int_b, b = slope_b, col = "blue")  # Lower bound
    abline(v = x_max, col = "green")  # x.max
    dev.off()
    
    # Print fit statistics
    message("\nData-driven bound statistics:")
    message(sprintf("Base fit: y = %.3fx + %.3f", base_slope, base_intercept))
    message(sprintf("Residual SD: %.3f", resid_sd))
    message(sprintf("R-squared: %.3f", summary(lm_fit)$r.squared))
    
    return(list(
        slope.r = slope_r,
        int.r = int_r,
        slope.b = slope_b,
        int.b = int_b,
        x.max = x_max,
        fit = lm_fit
    ))
}

# Calculate data-driven parameters
bounds <- calculate_data_driven_bounds(dist_matrix, nn_value)

# Create parameter combinations around the data-driven values
parameter_combinations <- list(
    # Base parameters from data: Uses exactly 2.5 SD from regression line
    # This is our default choice as it balances sensitivity and specificity
    list(slope.r = bounds$slope.r, 
         int.r = bounds$int.r, 
         slope.b = bounds$slope.b, 
         int.b = bounds$int.b,
         description = "Base parameters (2.5 SD from regression)"),
    
    # Slightly more stringent: Tightens bounds by 5%
    # Use when you need more conservative outlier detection
    list(slope.r = bounds$slope.r * 0.95, 
         int.r = bounds$int.r * 0.95,
         slope.b = bounds$slope.b * 1.05, 
         int.b = bounds$int.b * 1.05,
         description = "Stringent bounds (bounds tightened by 5%)"),
    
    # Slightly more lenient: Relaxes bounds by 5%
    # Use when you want to identify only the most extreme outliers
    list(slope.r = bounds$slope.r * 1.05, 
         int.r = bounds$int.r * 1.05,
         slope.b = bounds$slope.b * 0.95, 
         int.b = bounds$int.b * 0.95,
         description = "Lenient bounds (bounds relaxed by 5%)")
)

outliers_results <- list()
for(i in seq_along(parameter_combinations)) {
    params <- parameter_combinations[[i]]
    
    plot_file <- file.path("test", "test_results", "plots", "dimensionality_reduction", "outliers", sprintf("outliers_set%d.png", i))
    png(plot_file, width = 800, height = 600, res = 100)
    outliers <- knnOutliers(urd_object, 
                          nn.1 = 1,
                          nn.2 = round(nn_value/3),
                          x.max = bounds$x.max,
                          slope.r = params$slope.r,
                          int.r = params$int.r,
                          slope.b = params$slope.b,
                          int.b = params$int.b,
                          title = sprintf("%s\n(Parameter Set %d)", params$description, i))
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

# Compare results and select best parameters
# We prefer the base parameter set (Set 1) if it gives reasonable results (1-10% outliers)
best_result <- NULL
if(outliers_results[[1]]$percent >= 1 && outliers_results[[1]]$percent <= 10) {
    best_result <- outliers_results[[1]]
    message("\nSelected Parameter Set 1 (Base parameters) as it gives reasonable outlier percentage.")
} else {
    # If base parameters don't work, try others
    for(i in seq_along(outliers_results)) {
        if(outliers_results[[i]]$percent >= 1 && outliers_results[[i]]$percent <= 10) {
            best_result <- outliers_results[[i]]
            message(sprintf("\nSelected Parameter Set %d (%s) as it gives reasonable outlier percentage.", 
                          i, outliers_results[[i]]$description))
            break
        }
    }
}

# If no good result found, use most lenient set
if(is.null(best_result)) {
    best_result <- outliers_results[[length(outliers_results)]]
    message(sprintf("\nNo parameter set gave 1-10%% outliers. Using most lenient set (%s).", 
                   best_result$description))
}

best_outliers <- best_result$outliers

# Save outlier information with parameter set details
writeLines(c(
    sprintf("# Outlier Detection Parameters Used:"),
    sprintf("# %s", best_result$description),
    sprintf("# Slope (red): %.3f, Intercept (red): %.3f", best_result$params$slope.r, best_result$params$int.r),
    sprintf("# Slope (blue): %.3f, Intercept (blue): %.3f", best_result$params$slope.b, best_result$params$int.b),
    sprintf("# Outlier percentage: %.2f%%", best_result$percent),
    sprintf("# Total cells analyzed: %d", n_cells),
    sprintf("# Number of outlier cells: %d", length(best_outliers)),
    "",
    "# Dataset Statistics:",
    sprintf("# Minimum distance: %.2f", min(all_distances)),
    sprintf("# Median distance: %.2f", median(all_distances)),
    sprintf("# 95th percentile (x.max): %.2f", x_max_value),
    sprintf("# Maximum distance: %.2f", max(all_distances)),
    "",
    "# Cell IDs identified as outliers:",
    best_outliers
), "./test/test_results/outlier_detection_parameters.txt")

# Also save just the outlier cell IDs in a separate file
write.table(best_outliers, 
           file = "./test/test_results/outlier_cells.txt", 
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
saveRDS(urd_object, "./test/test_data/test_urd_object_with_dimred.rds")
saveRDS(urd_object_clean, "./test/test_data/test_urd_object_clean.rds")

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
          "./test/test_results/dimensionality_reduction/parameters.csv", 
          row.names = FALSE, 
          quote = FALSE)

message("\nDimensionality reduction analysis complete!")
message("Results saved to:")
message("- URD object with dimensionality reduction: ./test/test_data/test_urd_object_with_dimred.rds")
message("- Clean URD object: ./test/test_data/test_urd_object_clean.rds")
message("- Parameter summary: ./test/test_results/dimensionality_reduction/parameters.csv")
message("- Plots: ./test/test_results/plots/dimensionality_reduction/") 