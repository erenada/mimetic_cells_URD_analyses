# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(RColorBrewer)
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

# Find the most recent test URD object
test_files <- list.files("test_data", pattern = "test_urd_object_.*_with_var_genes\\.rds$", full.names = TRUE)
if (length(test_files) == 0) {
  stop("No test URD object with variable genes found. Please run run_test_linage_analysis.R first.")
}
latest_test_file <- test_files[which.max(file.info(test_files)$mtime)]

# Load the test URD object
message(sprintf("Loading test URD object from: %s", latest_test_file))
urd_object <- readRDS(latest_test_file)

# Create output directories
dir.create("test_results/plots/dimensionality_reduction", recursive = TRUE, showWarnings = FALSE)

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
# Try different mp.factor values (adjusted for smaller dataset)
mp_factors <- c(1.5, 2)
pca_results <- list()

safe_pdf("test_results/plots/dimensionality_reduction/pca_parameter_selection.pdf", {
  for(mp in mp_factors) {
    urd_object <- calcPCA(urd_object, mp.factor = mp)
    pca_results[[as.character(mp)]] <- length(urd_object@pca.sig)
    pcSDPlot(urd_object)
    title(main = sprintf("mp.factor = %.1f", mp))
  }
})

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

safe_pdf("test_results/plots/dimensionality_reduction/tsne_parameter_selection.pdf", {
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
})

# Store the last tSNE result (highest perplexity) for downstream analysis
set.seed(123)
urd_object <- calcTsne(urd_object, perplexity = max(perplexity_values), theta = 0.5)

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

set.seed(123)
urd_object <- graphClustering(urd_object, 
                           dim.use = "pca", 
                           num.nn = nn_values,
                           do.jaccard = TRUE, 
                           method = "Louvain")

# Plot clustering results
safe_pdf("test_results/plots/dimensionality_reduction/clustering_parameter_selection.pdf", {
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
})

# 4. kNN and Outlier Detection
# -----------------------------
message("\n4. kNN and Outlier Detection")
# Calculate optimal kNN value (adjusted for smaller dataset)
nn_value <- as.integer(ceiling(sqrt(n_cells)/2))  # Reduced from original calculation

message(sprintf("\nChosen parameters for kNN:"))
message(sprintf("nn_value: %d (based on dataset size)", nn_value))
message(sprintf("nn.2: %d (1/5 of nn_value)", as.integer(round(nn_value/5))))

# Try different nn values
nn_values <- as.integer(unique(round(c(nn_value * 0.5, nn_value, nn_value * 1.5))))

safe_pdf("test_results/plots/dimensionality_reduction/knn_parameter_selection.pdf", {
  outliers_per_nn <- list()
  for(nn in nn_values) {
    urd_object <- calcKNN(urd_object, nn = nn)
    outliers <- knnOutliers(urd_object, 
                          nn.1 = 1,
                          nn.2 = as.integer(round(nn/5)),
                          x.max = 40,        # Original value
                          slope.r = 1.1,     # Original value
                          int.r = 2.9,       # Original value
                          slope.b = 0.85,    # Original value
                          int.b = 10,        # Original value
                          title = sprintf("Outliers (nn=%d)", nn))
    outliers_per_nn[[as.character(nn)]] <- outliers
  }
})

# Compare results
message("\nOutlier detection results with different nn values:")
for(nn in names(outliers_per_nn)) {
    n_outliers <- length(outliers_per_nn[[nn]])
    percent_outliers <- 100 * n_outliers / n_cells
    message(sprintf("nn=%s: %d outliers (%.1f%%)", nn, n_outliers, percent_outliers))
}

# Use the chosen nn_value for final outlier detection
urd_object <- calcKNN(urd_object, nn = nn_value)
outliers <- knnOutliers(urd_object, 
                      nn.1 = 1,
                      nn.2 = as.integer(round(nn_value/5)),
                      x.max = 40,        # Original value
                      slope.r = 1.1,     # Original value
                      int.r = 2.9,       # Original value
                      slope.b = 0.85,    # Original value
                      int.b = 10,        # Original value
                      title = "Final Outlier Detection")

# Save outlier information
write.table(outliers, 
           file = "test_results/outlier_cells.txt", 
           row.names = FALSE, 
           col.names = FALSE, 
           quote = FALSE)

message(sprintf("\nIdentified %d outliers (%.1f%% of cells)", 
                length(outliers), 100 * length(outliers) / n_cells))

# Create clean subset
cells.keep <- setdiff(colnames(urd_object@logupx.data), outliers)
urd_object_clean <- urdSubset(urd_object, cells.keep = cells.keep)

# Save results
saveRDS(urd_object, "test_data/test_urd_object_with_dimred.rds")
saveRDS(urd_object_clean, "test_data/test_urd_object_clean.rds")

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
             sprintf("%.1f%%", 100 * length(outliers) / n_cells))
)
write.csv(parameter_summary, 
          "test_results/parameter_summary.csv", 
          row.names = FALSE, 
          quote = FALSE)

message("\nTest analysis complete!")
message("Parameter summary saved to 'test_results/parameter_summary.csv'") 