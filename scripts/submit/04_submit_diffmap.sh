#!/bin/bash

#SBATCH -J URD_diffmap
#SBATCH -p medium
#SBATCH -t 3-00:00:00
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -o ../../logs/diffmap_%j.out
#SBATCH -e ../../logs/diffmap_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eren_ada@hms.harvard.edu

# Create necessary directories
cd ../..
mkdir -p logs data results

# Load required modules
module load R/4.2.2

# Check if dimensionality reduction results exist
if [ ! -f "data/urd_object_clean.rds" ]; then
    echo "Error: Dimensionality reduction results not found. Please run dimensionality reduction first."
    exit 1
fi

# Run diffusion map analysis
echo "Starting diffusion map analysis..."
Rscript scripts/R/run_diffusion_map.R

if [ $? -eq 0 ]; then
    echo "Diffusion map analysis completed successfully."
else
    echo "Error: Diffusion map analysis failed."
    exit 1
fi 