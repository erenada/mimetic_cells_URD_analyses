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

# Create output directories
dir.create("test_results/plots/diffusion_map", recursive = TRUE, showWarnings = FALSE)

# Load the test URD object with dimensionality reduction results
message("Loading test URD object...")
urd_object <- readRDS("test_data/test_urd_object_with_dimred.rds")

# Parameter exploration for diffusion map
message("\nExploring different parameter combinations for diffusion map...")

# Different parameter combinations to try (adjusted for smaller dataset)
knn_values <- c(50, 75, 100)  # Reduced from original values
sigma_values <- c(4, 6, 8)    # Reduced range

# Store results for comparison
results <- list()

# Create plots for parameter exploration
safe_pdf("test_results/plots/diffusion_map/parameter_exploration.pdf", {
  # Try different parameter combinations
  for(knn in knn_values) {
    for(sigma in sigma_values) {
      message(sprintf("\nTrying knn=%d, sigma=%d", knn, sigma))
      
      # Calculate diffusion map with current parameters
      urd_object <- calcDM(urd_object,
                         knn = knn,
                         sigma.use = sigma,
                         distance = "euclidean",
                         verbose = TRUE)
      
      # Plot first 12 diffusion components (reduced from 18 for test)
      plotDimArray(object = urd_object, 
                  reduction.use = "dm", 
                  dims.to.plot = 1:12,
                  label = "stage",
                  plot.title = sprintf("knn=%d, sigma=%d", knn, sigma),
                  alpha = 0.3)
      
      # Store the diffusion map
      results[[sprintf("knn%d_sigma%d", knn, sigma)]] <- urd_object@dm
    }
  }
})

# Check if parameters were provided as command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 2) {
    chosen_knn <- as.numeric(args[1])
    chosen_sigma <- as.numeric(args[2])
    
    message(sprintf("\nUsing chosen parameters: knn=%d, sigma=%d", chosen_knn, chosen_sigma))
    
    # Calculate final diffusion map with chosen parameters
    urd_object <- calcDM(urd_object,
                        knn = chosen_knn,
                        sigma.use = chosen_sigma,
                        distance = "euclidean",
                        verbose = TRUE)
    
    # Plot final diffusion map components
    safe_pdf("test_results/plots/diffusion_map/final_diffusion_components.pdf", {
      plotDimArray(urd_object, 
                  reduction.use = "dm", 
                  dims.to.plot = 1:12,    # Reduced from 18 for test
                  label = "stage",
                  plot.title = sprintf("Final Diffusion Map (knn=%d, sigma=%d)", chosen_knn, chosen_sigma),
                  alpha = 0.3)
    })
    
    # Save the object with diffusion map
    saveRDS(urd_object, sprintf("test_data/test_urd_object_with_dm_knn%d_sigma%d.rds", chosen_knn, chosen_sigma))
    
    message("\nTest diffusion map analysis completed successfully!")
    message("\nNext steps:")
    message("1. Define root cells")
    message("2. Calculate pseudotime using floodPseudotime()")
    message("3. Process transition probabilities")
} else {
    message("\nNo parameters provided. Please examine the parameter exploration plots")
    message("and run the script again with your chosen parameters.")
    message("Example: Rscript run_test_diffusion_map.R 75 6")
} 