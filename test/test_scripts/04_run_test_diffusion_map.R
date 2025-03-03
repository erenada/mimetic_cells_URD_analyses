# Function to evaluate diffusion map quality
evaluate_dm_quality <- function(urd_object) {
    # Get the diffusion components
    dc_coords <- urd_object@dm@eigenvectors
    
    # Calculate metrics for each DC
    n_dcs <- min(18, ncol(dc_coords))
    dc_metrics <- matrix(0, nrow=3, ncol=n_dcs)
    rownames(dc_metrics) <- c("variance_explained", "entropy", "stage_separation")
    
    for(i in 1:n_dcs) {
        dc <- dc_coords[,i]
        
        # Variance explained - scale up by 100 to make it more meaningful
        dc_metrics[1,i] <- (var(dc) / sum(apply(dc_coords, 2, var))) * 100
        
        # Entropy (using fewer breaks for stability)
        breaks <- seq(min(dc), max(dc), length.out=20)
        hist_counts <- hist(dc, breaks=breaks, plot=FALSE)$counts
        props <- hist_counts / sum(hist_counts)
        # Normalize entropy to 0-1 range by dividing by log(number of bins)
        dc_metrics[2,i] <- -sum(props * log(props + 1e-10)) / log(length(props))
        
        # Stage separation - already normalized by sd(dc)
        stage_vars <- tapply(dc, urd_object@group.ids$stage, var)
        dc_metrics[3,i] <- mean(sqrt(stage_vars), na.rm=TRUE) / sd(dc)
    }
    
    # Calculate mean metrics
    mean_var_explained <- mean(dc_metrics[1,])
    mean_entropy <- mean(dc_metrics[2,])
    mean_stage_sep <- mean(dc_metrics[3,])
    
    # Combine metrics using weighted sum instead of multiplication
    # Weights prioritize variance explained and stage separation
    quality_score <- (0.4 * mean_var_explained) + 
                    (0.2 * mean_entropy) + 
                    (0.4 * mean_stage_sep)
    
    return(list(
        overall_score = quality_score,
        variance_explained = dc_metrics[1,],
        entropy = dc_metrics[2,],
        stage_separation = dc_metrics[3,]
    ))
}

# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(RColorBrewer)
})

# Create all necessary directories
dir.create("./test/test_data", recursive = TRUE, showWarnings = FALSE)
dir.create("./test/test_results/diffusion_map", recursive = TRUE, showWarnings = FALSE)
dir.create("./test/test_results/plots/diffusion_map", recursive = TRUE, showWarnings = FALSE)

# Find the most recent test URD object with clean data
test_files <- list.files("./test/test_data", pattern = "test_urd_object_clean\\.rds$", full.names = TRUE)
if (length(test_files) == 0) {
  stop("No clean test URD object found. Please run run_test_dimensionality_reduction.R first.")
}
latest_test_file <- test_files[which.max(file.info(test_files)$mtime)]

# Load the test URD object
message(sprintf("Loading clean URD object from: %s", latest_test_file))
urd_object <- readRDS(latest_test_file)

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
# Use sqrt(n_cells) as recommended by URD documentation, but adjusted for test data
min_knn <- 5  # Minimum knn to ensure stability
base_knn <- max(min_knn, ceiling(sqrt(n_cells)/2))  # Smaller for test data, but not too small
message(sprintf("\nKNN calculation:"))
message(sprintf("Base knn (max of %d and sqrt of cells / 2): %d", min_knn, base_knn))

# Try both auto-determined and calculated knn values
knn_values <- list(
    auto = max(min_knn, ceiling(sqrt(n_cells)/3)),  # Slightly larger than /4 for better stability
    calc = base_knn  # Use adjusted sqrt(n_cells) for test data
)

# Create lists to store results for each knn value
dm_results_by_knn <- list()
dm_metrics_by_knn <- list()

# Function to estimate optimal sigma range
estimate_sigma_range <- function(urd_object, knn_value) {
    # Calculate global distance statistics regardless of KNN mode
    # This helps us set reasonable sigma ranges even for auto KNN
    if(is.null(urd_object@knn)) {
        # Use a reasonable KNN for distance calculation even if we'll use auto KNN later
        temp_knn <- ceiling(sqrt(ncol(urd_object@logupx.data))/2)  # Smaller for test data
        urd_object <- calcKNN(urd_object, k = temp_knn)
    }
    distances <- urd_object@knn$nn.dists
    
    # Calculate global distance statistics
    global_median_dist <- median(distances)
    global_dist_sd <- sd(distances)
    
    # For numeric sigma values, calculate based on distances
    if(!is.null(knn_value)) {
        # Use the specific KNN distances
        distances_knn <- distances[,1:min(knn_value, ncol(distances))]
        median_dist <- median(distances_knn)
        dist_sd <- sd(distances_knn)
        
        # Define range based on distance distribution
        # More conservative range for test data
        sigma_min <- max(2.0, round(median_dist * 0.15, 2))   # Start from 15% of median distance
        sigma_max <- min(12.0, round(median_dist * 1.0, 2))   # Up to 100% of median distance
        
        # Generate sequence of numeric sigma values
        n_steps <- 3  # Fewer steps for test data since we're adding NULL and "local"
        sigma_values <- round(exp(seq(log(sigma_min), log(sigma_max), length.out = n_steps)), 2)
        
        # Add NULL and "local" options at the beginning
        sigma_values <- c(NULL, "local", sigma_values)
    } else {
        # For auto knn, use NULL, "local", and data-driven fixed values
        # Generate fixed values based on global distance statistics
        # More conservative range for test data
        sigma_min <- max(2.0, round(global_median_dist * 0.15, 2))  # 15% of median
        sigma_mid <- round(global_median_dist * 0.5, 2)             # 50% of median
        sigma_max <- min(12.0, round(global_median_dist * 1.0, 2))  # 100% of median
        
        # Create numeric values first
        numeric_values <- c(sigma_min, sigma_mid, sigma_max)
        numeric_values <- unique(numeric_values[!is.na(numeric_values)])  # Remove any duplicates or NAs
        
        # Add NULL and "local" options at the beginning
        sigma_values <- c(NULL, "local", numeric_values)
        
        median_dist <- global_median_dist
        dist_sd <- global_dist_sd
    }
    
    message("\nDistance statistics for sigma calculation:")
    if(!is.null(knn_value)) {
        message(sprintf("Using knn = %d", knn_value))
    } else {
        message("Using auto-determined knn (letting destiny optimize)")
        message(sprintf("Global distance statistics used for sigma range:"))
    }
    message(sprintf("Median distance: %.2f", median_dist))
    message(sprintf("Distance SD: %.2f", dist_sd))
    
    message("\nSigma values to try:")
    message("- NULL (destiny auto-detection)")
    message("- local (adaptive per cell)")
    
    # Get only the numeric values for percentage calculation
    numeric_sigmas <- as.numeric(sigma_values[!sapply(sigma_values, function(x) is.null(x) || is.character(x))])
    if(length(numeric_sigmas) > 0) {
        message(sprintf("- Fixed values: %s", paste(numeric_sigmas, collapse=", ")))
        message(sprintf("  (based on %.2f%% to %.2f%% of median distance)", 
               min(numeric_sigmas)/median_dist*100,
               max(numeric_sigmas)/median_dist*100))
    }
    
    return(list(
        sigma_values = sigma_values,
        median_dist = median_dist,
        dist_sd = dist_sd
    ))
}

# Try different knn values
for(knn_name in names(knn_values)) {
    knn_value <- knn_values[[knn_name]]
    message(sprintf("\nTrying %s knn value: %s", knn_name, if(is.null(knn_value)) "auto" else knn_value))
    
    # Calculate sigma range
    sigma_info <- estimate_sigma_range(urd_object, knn_value)
    sigma_values <- sigma_info$sigma_values
    
    message(sprintf("\nParameter selection for %s knn:", knn_name))
    message(sprintf("sigma values: %s", paste(sigma_values, collapse=", ")))
    
    # Create storage for this knn value
    dm_results_by_knn[[knn_name]] <- list()
    dm_metrics_by_knn[[knn_name]] <- list()
    
    # Calculate diffusion maps for each sigma value
    for(sigma in sigma_values) {
        # Convert sigma to string for messages
        sigma_str <- if(is.null(sigma)) "NULL" else if(is.character(sigma)) sigma else sprintf("%.2f", sigma)
        message(sprintf("\nCalculating diffusion map with knn = %s, sigma = %s...", 
                       if(is.null(knn_value)) "auto" else knn_value, sigma_str))
        
        # Calculate diffusion map with error handling
        dm_success <- tryCatch({
            urd_object <- calcDM(urd_object, 
                                knn = knn_value, 
                                sigma.use = sigma,
                                density.norm = TRUE,
                                verbose = TRUE)
            
            if (!is.null(urd_object@dm)) {
                TRUE
            } else {
                message("Diffusion map calculation returned NULL")
                FALSE
            }
        }, error = function(e) {
            message(sprintf("Error in diffusion map calculation: %s", e$message))
            FALSE
        })
        
        if(dm_success) {
            message(sprintf("Diffusion map calculation successful"))
            
            # Store results with combined key
            result_key <- sprintf("%s_%s", knn_name, 
                                if(is.null(sigma)) "auto" 
                                else if(is.character(sigma)) sigma 
                                else sprintf("%.2f", sigma))
            dm_results_by_knn[[knn_name]][[as.character(sigma)]] <- urd_object@dm
            
            # Evaluate quality
            tryCatch({
                quality_metrics <- evaluate_dm_quality(urd_object)
                dm_metrics_by_knn[[knn_name]][[as.character(sigma)]] <- quality_metrics
                
                # Save diffusion map with quality score in filename
                quality_score <- quality_metrics$overall_score
                detailed_filename <- sprintf(
                    "./test/test_results/diffusion_map/dm_%s_score%.3f.rds",
                    result_key, quality_score
                )

                # Save with basic metrics
                saveRDS(list(
                    diffusion_map = dm_results_by_knn[[knn_name]][[as.character(sigma)]],
                    metrics = list(
                        knn = knn_value,
                        sigma = as.numeric(sigma),
                        quality_score = quality_score
                    )
                ), detailed_filename)

                # Update summary dataframe
                if(!exists("all_metrics_summary")) {
                    all_metrics_summary <- data.frame()
                }

                all_metrics_summary <- rbind(all_metrics_summary,
                    data.frame(
                        knn = if(is.null(knn_value)) "auto" else knn_value,
                        sigma = as.numeric(sigma),
                        quality_score = quality_score,
                        filename = basename(detailed_filename)
                    )
                )
                
                # Create visualization
                png(sprintf("./test/test_results/plots/diffusion_map/dm_%s.png", result_key),
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
                            outer.title = sprintf("STAGE - Diffusion Map knn %s, Sigma %s", 
                                                knn_name,
                                                if(is.null(sigma)) "auto" 
                                                else if(is.character(sigma)) sigma 
                                                else sprintf("%.2f", sigma)),
                            legend = FALSE, 
                            alpha = 0.3, 
                            discrete.colors = stage_colors)
                
                dev.off()
                
                # Force garbage collection
                gc()
                
            }, error = function(e) {
                message(sprintf("Error in evaluating quality for %s knn, sigma %s: %s", 
                              knn_name, 
                              if(is.null(sigma)) "auto" 
                              else if(is.character(sigma)) sigma 
                              else sprintf("%.2f", sigma), 
                              e$message))
            })
        } else {
            message(sprintf("Diffusion map calculation failed for %s knn, sigma %s", 
                          knn_name,
                          if(is.null(sigma)) "auto" 
                          else if(is.character(sigma)) sigma 
                          else sprintf("%.2f", sigma)))
        }
    }
}

# Report final status and quality metrics
message("\nDiffusion Map Analysis Summary:")
message(sprintf("Attempted knn values: %s", paste(names(knn_values), collapse=", ")))
message(sprintf("Successful calculations: %s", 
                paste(names(dm_results_by_knn), collapse=", ")))

# Find best combination across all knn values
best_quality <- -Inf
best_knn <- NULL
best_sigma <- NULL

if(length(dm_metrics_by_knn) > 0) {
    message("\nQuality metrics for each knn value:")
    for(knn_name in names(dm_metrics_by_knn)) {
        message(sprintf("\nknn: %s", knn_name))
        for(sigma in names(dm_metrics_by_knn[[knn_name]])) {
            metrics <- dm_metrics_by_knn[[knn_name]][[sigma]]
            quality_score <- metrics$overall_score
            # Format sigma string for display
            sigma_str <- if(is.null(sigma)) "auto" 
                        else if(is.character(sigma)) sigma 
                        else sprintf("%.2f", as.numeric(sigma))
            message(sprintf("\nSigma: %s", sigma_str))
            message(sprintf("  Overall quality score: %.3f", quality_score))
            message(sprintf("  Mean variance explained: %.3f", mean(metrics$variance_explained)))
            message(sprintf("  Mean entropy: %.3f", mean(metrics$entropy)))
            message(sprintf("  Mean stage separation: %.3f", mean(metrics$stage_separation)))
            
            # Update best if better quality found
            if(quality_score > best_quality) {
                best_quality <- quality_score
                best_knn <- knn_name
                best_sigma <- sigma  # Store original sigma value, not converted
            }
        }
    }
    
    # Format best sigma for display
    best_sigma_str <- if(is.null(best_sigma)) "auto"
                      else if(is.character(best_sigma)) best_sigma
                      else sprintf("%.2f", as.numeric(best_sigma))
    
    message(sprintf("\nBest combination found:"))
    message(sprintf("KNN: %s", best_knn))
    message(sprintf("Sigma: %s", best_sigma_str))
    message(sprintf("Quality score: %.3f", best_quality))
    
    # Use the best combination - use original sigma value
    urd_object@dm <- dm_results_by_knn[[best_knn]][[as.character(best_sigma)]]
} else {
    stop("No successful diffusion map calculations. Cannot proceed.")
}

# Save the URD object with the final diffusion map
saveRDS(urd_object, "./test/test_data/test_urd_object_with_dm.rds")

# Save parameter selection summary
parameter_summary <- data.frame(
    parameter = c("Dataset size", 
                 "Number of stages",
                 "Complexity factor",
                 "Best knn",
                 "Best sigma",
                 "Quality score"),
    value = c(sprintf("%d cells", n_cells),
              sprintf("%d", n_stages),
              sprintf("%.2f", complexity_factor),
              best_knn,
              best_sigma_str,
              sprintf("%.3f", best_quality))
)
write.csv(parameter_summary, 
          "./test/test_results/diffusion_map/parameters.csv", 
          row.names = FALSE, 
          quote = FALSE)

# Save the summary of all diffusion maps
write.csv(all_metrics_summary, 
          "./test/test_results/diffusion_map/all_diffusion_maps_summary.csv",
          row.names = FALSE)

message("\nDiffusion map analysis complete!")
message("Results saved to:")
message("- URD object: ./test/test_data/test_urd_object_with_dm.rds")
message("- Parameter summary: ./test/test_results/diffusion_map/parameters.csv")
message("- Individual maps: ./test/test_results/diffusion_map/")
message("- Plots: ./test/test_results/plots/diffusion_map/")