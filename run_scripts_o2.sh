#!/bin/bash
#SBATCH -c 20                               # Request 20 cores
#SBATCH -t 4-00:00:00                       # Runtime in D-HH:MM format (4 days)
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=64G                           # Memory total in MiB (for all cores)
#SBATCH -o logs/urd_%j.out                  # File to which STDOUT will be written, %j is job ID
#SBATCH -e logs/urd_%j.err                  # File to which STDERR will be written, %j is job ID
#SBATCH --mail-user=eren_ada@hms.harvard.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=URD_dimred              # Job name updated to reflect actual task

# Create necessary directories
mkdir -p logs results/plots/dimensionality_reduction

# Change to the working directory
cd /n/groups/immdiv-bioinfo/eren/mimetic_cells_URD_analyses

# Load required modules
module load gcc/9.2.0
module load R/4.4.0

# Print some information about the job
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Starting time: $(date)"
echo "Working directory: $(pwd)"

# Verify URD object exists
if [ ! -f "data/initial_urd_object_20250210_1206.rds" ]; then
    echo "Error: URD object file not found at data/initial_urd_object_20250210_1206.rds"
    exit 1
fi

# Run dimensionality reduction analysis
echo "Starting dimensionality reduction analysis..."

echo "Running dimensionality reduction..."
Rscript run_dimensionality_reduction.R
if [ $? -ne 0 ]; then
    echo "Error in run_dimensionality_reduction.R"
    exit 1
fi

echo "Analysis completed at: $(date)"