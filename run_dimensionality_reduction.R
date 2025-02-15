# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(RColorBrewer)
})
dev.off()
# Load the URD object with variable genes
urd_object <- readRDS("data/initial_urd_object_20250210_1206.rds")

# Create output directories
dir.create("results/plots/dimensionality_reduction", recursive = TRUE, showWarnings = FALSE)

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
# Try different mp.factor values
mp_factors <- c(1.5, 2, 2.5)
pca_results <- list()

pdf("results/plots/dimensionality_reduction/pca_parameter_selection.pdf")
for(mp in mp_factors) {
  urd_object <- calcPCA(urd_object, mp.factor = mp)
  pca_results[[as.character(mp)]] <- length(urd_object@pca.sig)
  pcSDPlot(urd_object)
  title(main = sprintf("mp.factor = %.1f", mp))
}
dev.off()

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
# Rule of thumb: perplexity should be between 5 and sqrt(n_cells)/3
min_perp <- 5
max_perp <- min(50, floor(sqrt(n_cells)/3))
perplexity_values <- unique(round(c(
    min_perp,
    floor(max_perp/2),
    max_perp
)))

message(sprintf("Testing perplexity values: %s", paste(perplexity_values, collapse=", ")))

pdf("results/plots/dimensionality_reduction/tsne_parameter_selection.pdf")
for(perp in perplexity_values) {
  message(sprintf("Computing tSNE with perplexity %d...", perp))
  set.seed(123)
  urd_object <- calcTsne(urd_object, perplexity = perp, theta = 0.5)
  
  # Plot stages
  plotDim(urd_object, 
          "stage",
          reduction.use = "tsne",
          plot.title = sprintf("tSNE of Stages (perplexity=%d)", perp),
          legend = TRUE)
}
dev.off()

# Store the last tSNE result (highest perplexity) for downstream analysis
set.seed(123)
urd_object <- calcTsne(urd_object, perplexity = max(perplexity_values), theta = 0.5)

# 3. Clustering Parameters
# -----------------------------
message("\n3. Graph-based Clustering")

# Calculate dataset complexity metrics
n_stages <- length(unique(urd_object@group.ids$stage))
cells_per_stage <- table(urd_object@group.ids$stage)
min_cells_per_stage <- min(cells_per_stage)
complexity_factor <- n_stages / log2(n_cells)  # Higher for more complex datasets

message("\nDataset complexity metrics:")
message(sprintf("Number of stages: %d", n_stages))
message(sprintf("Minimum cells per stage: %d", min_cells_per_stage))
message(sprintf("Complexity factor: %.2f", complexity_factor))

# Calculate number of nearest neighbors based on dataset size AND complexity
base_min_nn <- ceiling(sqrt(n_cells)/2)
base_max_nn <- ceiling(sqrt(n_cells))

# Adjust nn based on complexity:
# - For more complex data (many stages relative to cells), use fewer neighbors
# - For simpler data (few stages relative to cells), use more neighbors
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

set.seed(123)
urd_object <- graphClustering(urd_object, 
                           dim.use = "pca", 
                           num.nn = nn_values,
                           do.jaccard = TRUE, 
                           method = "Louvain")

# Plot clustering results
pdf("results/plots/dimensionality_reduction/clustering_parameter_selection.pdf")
for(nn in nn_values) {
  # Plot on PCA
  plotDim(urd_object, 
          sprintf("Louvain-%d", nn), 
          reduction.use = "pca",
          legend = TRUE,
          plot.title = sprintf("Louvain-Jaccard Clustering on PCA (%d NNs)", nn))
  
  # Plot on tSNE
  plotDim(urd_object, 
          sprintf("Louvain-%d", nn), 
          reduction.use = "tsne",
          legend = TRUE,
          plot.title = sprintf("Louvain-Jaccard Clustering on tSNE (%d NNs)", nn))
}
dev.off()

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

message(sprintf("\nChosen parameters for kNN:"))
message(sprintf("nn_value: %d (based on dataset size)", nn_value))
message(sprintf("nn.2: %d (1/3 of nn_value)", round(nn_value/3)))

# Try different parameter combinations
parameter_combinations <- list(
    list(slope.r = 1.1, int.r = 2.9, slope.b = 0.85, int.b = 10),  # Original
    list(slope.r = 1.2, int.r = 3.5, slope.b = 0.9, int.b = 8),    # More lenient
    list(slope.r = 1.3, int.r = 4.0, slope.b = 0.95, int.b = 6)    # Most lenient
)

pdf("results/plots/dimensionality_reduction/outlier_parameter_selection.pdf")
outliers_results <- list()
for(i in seq_along(parameter_combinations)) {
    params <- parameter_combinations[[i]]
    
    urd_object <- calcKNN(urd_object, nn = nn_value)
    outliers <- knnOutliers(urd_object, 
                          nn.1 = 1,
                          nn.2 = round(nn_value/3),
                          x.max = 40,
                          slope.r = params$slope.r,
                          int.r = params$int.r,
                          slope.b = params$slope.b,
                          int.b = params$int.b,
                          title = sprintf("Parameter Set %d", i))
    
    outliers_results[[i]] <- list(
        params = params,
        outliers = outliers,
        percent = 100 * length(outliers) / n_cells
    )
}
dev.off()

# Compare results
message("\nOutlier detection results with different parameter sets:")
for(i in seq_along(outliers_results)) {
    result <- outliers_results[[i]]
    message(sprintf("Parameter Set %d:", i))
    message(sprintf("  Parameters: slope.r=%.1f, int.r=%.1f, slope.b=%.2f, int.b=%.1f",
                   result$params$slope.r, result$params$int.r,
                   result$params$slope.b, result$params$int.b))
    message(sprintf("  Outliers: %d (%.1f%%)", 
                   length(result$outliers), result$percent))
}

# Use the parameter set that gives reasonable outlier percentage (between 1% and 10%)
best_params <- NULL
best_outliers <- NULL
for(result in outliers_results) {
    if(result$percent >= 1 && result$percent <= 10) {
        best_params <- result$params
        best_outliers <- result$outliers
        break
    }
}

if(is.null(best_params)) {
    # If no parameter set gives desired range, use most lenient set
    best_params <- parameter_combinations[[length(parameter_combinations)]]
    urd_object <- calcKNN(urd_object, nn = nn_value)
    best_outliers <- knnOutliers(urd_object, 
                               nn.1 = 1,
                               nn.2 = round(nn_value/3),
                               x.max = 40,
                               slope.r = best_params$slope.r,
                               int.r = best_params$int.r,
                               slope.b = best_params$slope.b,
                               int.b = best_params$int.b,
                               title = "Final Outlier Detection")
}

# Save outlier information
write.table(best_outliers, 
           file = "results/outlier_cells.txt", 
           row.names = FALSE, 
           col.names = FALSE, 
           quote = FALSE)

message(sprintf("\nIdentified %d outliers (%.1f%% of cells)", 
                length(best_outliers), 100 * length(best_outliers) / n_cells))

# Create clean subset
cells.keep <- setdiff(colnames(urd_object@logupx.data), best_outliers)
urd_object_clean <- urdSubset(urd_object, cells.keep = cells.keep)

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
          "results/parameter_summary.csv", 
          row.names = FALSE, 
          quote = FALSE)

message("\nAnalysis complete!")
message("Parameter summary saved to 'results/parameter_summary.csv'") 