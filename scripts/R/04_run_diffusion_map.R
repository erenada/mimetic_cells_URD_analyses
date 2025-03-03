# load libraries
suppressPackageStartupMessages({
  library(URD)
  library(Seurat)
  library(RColorBrewer)
})

# Create all necessary directories
dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/diffusion_map", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots/diffusion_map", recursive = TRUE, showWarnings = FALSE)

# Find the most recent URD object with clean data
urd_object <- readRDS("data/urd_object_clean.rds")

# Create output directories
dir.create("results/plots/diffusion_map", recursive = TRUE, showWarnings = FALSE)

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
# Use sqrt(n_cells) as recommended by URD documentation
min_knn <- 10  # Minimum knn to ensure stability
base_knn <- max(min_knn, ceiling(sqrt(n_cells)))  # Use full sqrt(n_cells) as recommended
message(sprintf("\nKNN calculation:"))
message(sprintf("Base knn (sqrt of cells): %d", base_knn))

# Try both auto-determined and calculated knn values
knn_values <- list(
    auto = NULL,  # Let destiny determine optimal knn as recommended in docs
    calc = base_knn  # Use sqrt(n_cells) as recommended
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
        temp_knn <- ceiling(sqrt(ncol(urd_object@logupx.data)))
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
        sigma_min <- max(2.0, round(median_dist * 0.3, 2))   # Start from 30% of median distance
        sigma_max <- min(20.0, round(median_dist * 1.5, 2))  # Up to 150% of median distance
        
        # Generate sequence of numeric sigma values
        n_steps <- 4  # Reduced steps since we're adding NULL and "local"
        sigma_values <- round(exp(seq(log(sigma_min), log(sigma_max), length.out = n_steps)), 2)
        
        # Add NULL and "local" options at the beginning
        sigma_values <- c(NULL, "local", sigma_values)
    } else {
        # For auto knn, use NULL, "local", and data-driven fixed values
        # Generate fixed values based on global distance statistics
        sigma_min <- max(2.0, round(global_median_dist * 0.3, 2))  # 30% of median
        sigma_mid1 <- round(global_median_dist * 0.75, 2)          # 75% of median
        sigma_mid2 <- round(global_median_dist * 1.0, 2)           # 100% of median
        sigma_max <- min(20.0, round(global_median_dist * 1.5, 2)) # 150% of median
        
        # Create numeric values first
        numeric_values <- c(sigma_min, sigma_mid1, sigma_mid2, sigma_max)
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
                dm_metrics_by_knn[[knn_name]][[as.character(sigma)]] <- evaluate_dm_quality(urd_object)
                
                # Save diffusion map with quality score in filename
                quality_score <- dm_metrics_by_knn[[knn_name]][[as.character(sigma)]]$overall_score
                detailed_filename <- sprintf(
                    "results/diffusion_map/dm_%s_score%.3f.rds",
                    result_key, quality_score
                )

                # Save with basic metrics - use result_key for accessing dm_results_by_knn
                saveRDS(list(
                    diffusion_map = dm_results_by_knn[[knn_name]][[as.character(sigma)]],
                    metrics = list(
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
                        sigma = as.numeric(sigma),
                        quality_score = quality_score,
                        filename = basename(detailed_filename)
                    )
                )
                
                # Create visualization
                png(sprintf("results/plots/diffusion_map/dm_%s.png", result_key),
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
            
            # Update best if better quality found
            if(quality_score > best_quality) {
                best_quality <- quality_score
                best_knn <- knn_name
                best_sigma <- sigma  # Store original sigma value, not converted
            }
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
            # Format sigma string for display
            sigma_str <- if(is.null(sigma)) "auto" 
                        else if(is.character(sigma)) sigma 
                        else sprintf("%.2f", as.numeric(sigma))
            message(sprintf("\nSigma: %s", sigma_str))
            message(sprintf("  Overall quality score: %.3f", metrics$overall_score))
            message(sprintf("  Mean variance explained: %.3f", mean(metrics$variance_explained)))
            message(sprintf("  Mean entropy: %.3f", mean(metrics$entropy)))
            message(sprintf("  Mean stage separation: %.3f", mean(metrics$stage_separation)))
            
            # Update best if better quality found
            if(metrics$overall_score > best_quality) {
                best_quality <- metrics$overall_score
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
saveRDS(urd_object, "data/urd_object_with_dm.rds")

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
              sprintf("%s", names(knn_values)),
              sprintf("%.2f-%.2f", min(sapply(knn_values, function(x) if(is.null(x)) NA else x)),
              sprintf("%.2f", best_sigma),
              sprintf("%.3f", best_quality))
)
write.csv(parameter_summary, 
          "results/diffusion_map/parameters.csv", 
          row.names = FALSE, 
          quote = FALSE)

message("\nDiffusion map analysis complete!")
message("Results saved to:")
message("- URD object: data/urd_object_with_dm.rds")
message("- Parameter summary: results/diffusion_map/parameters.csv")
message("- Individual maps: results/diffusion_map/")
message("- Plots: results/plots/diffusion_map/") 