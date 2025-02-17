#!/bin/bash

#SBATCH -J URD_create
#SBATCH -p short
#SBATCH -t 0-12:00:00
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -o ../../logs/create_urd_%j.out
#SBATCH -e ../../logs/create_urd_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eren_ada@hms.harvard.edu

# Create necessary directories
cd ../..
mkdir -p logs data results

# Load required modules
module load R/4.2.2

# Run URD object creation
echo "Starting URD object creation..."
Rscript scripts/R/01_create_urd_object.R

if [ $? -eq 0 ]; then
    echo "URD object creation completed successfully."
else
    echo "Error: URD object creation failed."
    exit 1
fi 