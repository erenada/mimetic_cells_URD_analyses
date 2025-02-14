#!/bin/bash

# Create necessary directories
mkdir -p test_data test_results/plots logs

# Set up logging
exec 1> >(tee "logs/test_pipeline_$(date +%Y%m%d_%H%M%S).log")
exec 2>&1

echo "Starting test pipeline at $(date)"

# Step 1: Create test URD object
echo -e "\n=== Step 1: Creating test URD object ==="
Rscript create_test_urd_object.R
if [ $? -ne 0 ]; then
    echo "Error in create_test_urd_object.R"
    exit 1
fi

# Step 2: Run dimensionality reduction
echo -e "\n=== Step 2: Running dimensionality reduction ==="
Rscript run_test_dimensionality_reduction.R
if [ $? -ne 0 ]; then
    echo "Error in run_test_dimensionality_reduction.R"
    exit 1
fi

# Step 3: Run diffusion map parameter exploration
echo -e "\n=== Step 3: Running diffusion map parameter exploration ==="
Rscript run_test_diffusion_map.R
if [ $? -ne 0 ]; then
    echo "Error in run_test_diffusion_map.R"
    exit 1
fi

echo -e "\nTest pipeline completed at $(date)"
echo "Please examine the results in test_results/plots/"
echo "Once satisfied with parameter exploration, run the final diffusion map step with chosen parameters:"
echo "Rscript run_test_diffusion_map.R <knn> <sigma>" 