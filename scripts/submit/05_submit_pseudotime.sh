#!/bin/bash

#SBATCH -J URD_pseudotime
#SBATCH -p medium
#SBATCH -t 3-00:00:00
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -o ../../logs/pseudotime_%j.out
#SBATCH -e ../../logs/pseudotime_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eren_ada@hms.harvard.edu

# Create necessary directories
cd ../..
mkdir -p logs data results

# Load required modules
module load R/4.2.2

# Check if diffusion map results exist
if [ ! -f "data/urd_object_with_dm.rds" ]; then
    echo "Error: Diffusion map results not found. Please run diffusion map analysis first."
    exit 1
fi

# Run pseudotime analysis
echo "Starting pseudotime analysis..."
Rscript scripts/R/run_pseudotime.R

if [ $? -eq 0 ]; then
    echo "Pseudotime analysis completed successfully."
else
    echo "Error: Pseudotime analysis failed."
    exit 1
fi 