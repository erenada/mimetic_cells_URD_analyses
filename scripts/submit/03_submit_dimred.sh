#!/bin/bash

#SBATCH -J URD_dimred
#SBATCH -p medium
#SBATCH -t 3-00:00:00
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -o ../../logs/dimred_%j.out
#SBATCH -e ../../logs/dimred_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eren_ada@hms.harvard.edu

# Create necessary directories
cd ../..
mkdir -p logs data results

# Load required modules
module load R/4.2.2

# Check if variable genes analysis results exist
if [ ! -f "data/urd_object_with_var_genes.rds" ]; then
    echo "Error: Variable genes analysis results not found. Please run find_variable_genes.R first."
    exit 1
fi

# Run dimensionality reduction
echo "Starting dimensionality reduction analysis..."
Rscript scripts/R/run_dimensionality_reduction.R

if [ $? -eq 0 ]; then
    echo "Dimensionality reduction completed successfully."
else
    echo "Error: Dimensionality reduction failed."
    exit 1
fi 