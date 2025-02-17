# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(RColorBrewer)
})

# Find the most recent test URD object with clean data
test_files <- list.files("test_data", pattern = "test_urd_object_clean\\.rds$", full.names = TRUE)
if (length(test_files) == 0) {
  stop("No clean test URD object found. Please run run_test_dimensionality_reduction.R first.")
}
latest_test_file <- test_files[which.max(file.info(test_files)$mtime)]

# Load the test URD object
message(sprintf("Loading clean URD object from: %s", latest_test_file))
urd_object <- readRDS(latest_test_file)

# Create output directories
dir.create("test_results/plots/diffusion_map", recursive = TRUE, showWarnings = FALSE)

message("Starting diffusion map calculation...")

# Get dataset characteristics for parameter selection
n_cells <- ncol(urd_object@logupx.data)
n_stages <- length(unique(urd_object@group.ids$stage))
cells_per_stage <- table(urd_object@group.ids$stage)
min_cells_per_stage <- min(cells_per_stage)
complexity_factor <- n_stages / log2(n_cells)

message(sprintf("\nDataset characteristics:"))
message(sprintf("Number of cells: %d", n_cells))
message(sprintf("Number of stages: %d", n_stages))
message(sprintf("Minimum cells per stage: %d", min_cells_per_stage))
message(sprintf("Complexity factor: %.2f", complexity_factor))

# Calculate optimal knn based on data characteristics
# Base calculation considering dataset size and complexity
base_knn <- ceiling(sqrt(n_cells)/2)
complexity_adjustment <- 1 / (1 + log2(complexity_factor))
knn_value <- max(50, min(200, ceiling(base_knn * complexity_adjustment)))

# Function to estimate optimal sigma range
estimate_sigma_range <- function(urd_object, knn_value) {
    # Calculate distances to nearest neighbors
    if(is.null(urd_object@knn)) {
        urd_object <- calcKNN(urd_object, k = knn_value)
    }
    distances <- urd_object@knn$nn.dists[,1:min(knn_value, ncol(urd_object@knn$nn.dists))]
    
    # Calculate statistics
    median_dist <- median(distances)
    dist_sd <- sd(distances)
    
    # Define range based on distance distribution
    # Lower bound: captures local structure
    # Upper bound: allows for broader connections without over-smoothing
    sigma_min <- max(2, round(median_dist - 1.5 * dist_sd, 1))
    sigma_max <- min(12, round(median_dist + 1.5 * dist_sd, 1))
    
    # Generate sequence of sigma values
    n_steps <- 4  # Number of sigma values to test
    sigma_values <- round(seq(sigma_min, sigma_max, length.out = n_steps), 1)
    
    return(list(
        sigma_values = sigma_values,
        median_dist = median_dist,
        dist_sd = dist_sd
    ))
}

# Calculate sigma range
sigma_info <- estimate_sigma_range(urd_object, knn_value)
sigma_values <- sigma_info$sigma_values

message(sprintf("\nParameter selection (data-driven approach):"))
message(sprintf("knn: %d (adjusted for dataset size and complexity)", knn_value))
message(sprintf("sigma values: %s", paste(sigma_values, collapse=", ")))
message(sprintf("  - Based on distance distribution:"))
message(sprintf("    Median distance: %.2f", sigma_info$median_dist))
message(sprintf("    Distance SD: %.2f", sigma_info$dist_sd))

# Create a list to store diffusion maps and their quality metrics
dm_results <- list()
dm_metrics <- list()

# Function to evaluate diffusion map quality
evaluate_dm_quality <- function(urd_object) {
    # Get the diffusion components
    dc_coords <- urd_object@dm@eigenvectors  # Fixed: correct way to access DCs
    
    # Calculate metrics for each DC
    n_dcs <- min(18, ncol(dc_coords))
    dc_metrics <- matrix(0, nrow=3, ncol=n_dcs)
    rownames(dc_metrics) <- c("variance_explained", "entropy", "stage_separation")
    
    for(i in 1:n_dcs) {
        dc <- dc_coords[,i]
        
        # Variance explained
        dc_metrics[1,i] <- var(dc) / sum(apply(dc_coords, 2, var))
        
        # Entropy (using fewer breaks for stability)
        breaks <- seq(min(dc), max(dc), length.out=20)
        hist_counts <- hist(dc, breaks=breaks, plot=FALSE)$counts
        props <- hist_counts / sum(hist_counts)
        dc_metrics[2,i] <- -sum(props * log(props + 1e-10))
        
        # Stage separation
        stage_vars <- tapply(dc, urd_object@group.ids$stage, var)
        dc_metrics[3,i] <- mean(sqrt(stage_vars), na.rm=TRUE) / sd(dc)
    }
    
    # Combine metrics
    quality_score <- mean(dc_metrics[1,]) * mean(dc_metrics[2,]) * mean(dc_metrics[3,])
    
    return(list(
        overall_score = quality_score,
        variance_explained = dc_metrics[1,],
        entropy = dc_metrics[2,],
        stage_separation = dc_metrics[3,]
    ))
}

# Calculate diffusion maps for each sigma value
for(sigma in sigma_values) {
    message(sprintf("\nCalculating diffusion map with sigma = %.1f...", sigma))
    
    # Calculate diffusion map with error handling
    dm_success <- tryCatch({
        urd_object <- calcDM(urd_object, knn = knn_value, sigma.use = sigma)
        !is.null(urd_object@dm)
    }, error = function(e) {
        message(sprintf("Error in diffusion map calculation for sigma %.1f: %s", sigma, e$message))
        FALSE
    })
    
    if(dm_success) {
        message(sprintf("Diffusion map calculation successful for sigma %.1f", sigma))
        
        # Store the diffusion map
        dm_results[[as.character(sigma)]] <- urd_object@dm
        
        # Evaluate quality
        dm_metrics[[as.character(sigma)]] <- evaluate_dm_quality(urd_object)
        
        # Save individual diffusion map
        saveRDS(urd_object@dm, 
                file = sprintf("test_results/dm_sigma%.1f.rds", sigma))
        
        # Create visualization
        png(sprintf("test_results/plots/diffusion_map/dm_sigma%.1f.png", sigma),
            width = 1200, height = 1200, res = 150)
        
        # Stage color palette (expanded for more stages)
        n_stages <- length(unique(urd_object@group.ids$stage))
        stage_colors <- c(
            "#CCCCCC",  # First color
            RColorBrewer::brewer.pal(9, "Set1"),  # 9 colors
            RColorBrewer::brewer.pal(12, "Paired"),  # 12 colors
            RColorBrewer::brewer.pal(8, "Dark2")  # 8 more colors
        )[1:n_stages]  # Take only as many colors as we have stages
        
        # Plot first 18 diffusion components
        plotDimArray(object = urd_object, 
                    reduction.use = "dm", 
                    dims.to.plot = 1:18, 
                    label = "stage", 
                    plot.title = "", 
                    outer.title = sprintf("STAGE - Diffusion Map Sigma %.1f", sigma),
                    legend = FALSE, 
                    alpha = 0.3, 
                    discrete.colors = stage_colors)
        
        dev.off()
        
        # Force garbage collection
        gc()
        
    } else {
        message(sprintf("Diffusion map calculation failed for sigma %.1f", sigma))
    }
}

# Report final status and quality metrics
message("\nDiffusion Map Analysis Summary:")
message(sprintf("Attempted sigma values: %s", paste(sigma_values, collapse=", ")))
message(sprintf("Successful calculations: %s", 
                paste(names(dm_results), collapse=", ")))

if(length(dm_metrics) > 0) {
    message("\nQuality metrics for each sigma value:")
    for(sigma in names(dm_metrics)) {
        metrics <- dm_metrics[[sigma]]
        message(sprintf("\nSigma %.1f:", as.numeric(sigma)))
        message(sprintf("  Overall quality score: %.3f", metrics$overall_score))
        message(sprintf("  Mean variance explained: %.3f", mean(metrics$variance_explained)))
        message(sprintf("  Mean entropy: %.3f", mean(metrics$entropy)))
        message(sprintf("  Mean stage separation: %.3f", mean(metrics$stage_separation)))
    }
    
    # Select optimal sigma based on quality metrics
    quality_scores <- sapply(dm_metrics, function(x) x$overall_score)
    optimal_sigma <- as.numeric(names(quality_scores)[which.max(quality_scores)])
    
    message(sprintf("\nSelected optimal sigma %.1f based on quality metrics", optimal_sigma))
    
    # Use the optimal sigma result
    urd_object@dm <- dm_results[[as.character(optimal_sigma)]]
} else {
    stop("No successful diffusion map calculations. Cannot proceed.")
}

# Save the URD object with the final diffusion map
saveRDS(urd_object, "test_data/test_urd_object_with_dm.rds")

# Save parameter selection summary
parameter_summary <- data.frame(
    parameter = c("Dataset size", 
                 "Number of stages",
                 "Complexity factor",
                 "Selected knn",
                 "Tested sigma range",
                 "Optimal sigma",
                 "Quality score"),
    value = c(sprintf("%d cells", n_cells),
              sprintf("%d", n_stages),
              sprintf("%.2f", complexity_factor),
              sprintf("%d", knn_value),
              sprintf("%.1f-%.1f", min(sigma_values), max(sigma_values)),
              sprintf("%.1f", optimal_sigma),
              sprintf("%.3f", max(quality_scores)))
)
write.csv(parameter_summary, 
          "test_results/diffusion_map_parameters.csv", 
          row.names = FALSE, 
          quote = FALSE)

message("\nDiffusion map analysis complete!")
message("Results saved to 'test_data/test_urd_object_with_dm.rds'")
message("Parameter summary saved to 'test_results/diffusion_map_parameters.csv'")
message("Plots saved in 'test_results/plots/diffusion_map/'") 