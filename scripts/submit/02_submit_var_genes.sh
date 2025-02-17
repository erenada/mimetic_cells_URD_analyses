#!/bin/bash

#SBATCH -J URD_var_genes
#SBATCH -p short
#SBATCH -t 0-12:00:00
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -o ../../logs/var_genes_%j.out
#SBATCH -e ../../logs/var_genes_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eren_ada@hms.harvard.edu

# Create necessary directories
cd ../..
mkdir -p logs data results

# Load required modules
module load R/4.2.2

# Run variable genes analysis
echo "Starting variable genes analysis..."
Rscript scripts/R/find_variable_genes.R

if [ $? -eq 0 ]; then
    echo "Variable genes analysis completed successfully."
else
    echo "Error: Variable genes analysis failed."
    exit 1
fi 