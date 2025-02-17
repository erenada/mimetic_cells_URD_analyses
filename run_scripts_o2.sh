#!/bin/bash
#SBATCH -c 20                               # Request 20 cores
#SBATCH -t 12:00:00                       # Runtime in D-HH:MM format (12 hours)
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=64G                           # Memory total in MiB (for all cores)
#SBATCH -o logs/urd_%j.out                  # File to which STDOUT will be written, %j is job ID
#SBATCH -e logs/urd_%j.err                  # File to which STDERR will be written, %j is job ID
#SBATCH --mail-user=eren_ada@hms.harvard.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=URD_dimred              # Job name updated to reflect actual task

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to get memory usage
get_memory_usage() {
    free -h | awk 'NR==2{printf "Memory Usage: %s/%s (%.2f%%)\n", $3,$2,$3*100/$2 }'
}

# Function to check file creation and size
check_file() {
    local file="$1"
    if [ -f "$file" ]; then
        FILE_SIZE=$(stat --format=%s "$file")
        FILE_SIZE_MB=$(echo "scale=2; $FILE_SIZE/1024/1024" | bc)
        log_message "✓ $file successfully created ($FILE_SIZE_MB MB)"
        return 0
    else
        log_message "✗ $file not found"
        return 1
    fi
}

# Create necessary directories
mkdir -p logs results/plots/dimensionality_reduction

# Change to the working directory
cd /n/groups/immdiv-bioinfo/eren/mimetic_cells_URD_analyses

# Record start time
START_TIME=$(date +%s)

# Log system information
log_message "=== Job Information ==="
log_message "Job ID: $SLURM_JOB_ID"
log_message "Node: $SLURM_JOB_NODELIST"
log_message "Number of cores: $SLURM_CPUS_PER_TASK"
log_message "Working directory: $(pwd)"

# Log system specifications
log_message "=== System Information ==="
log_message "CPU Info:"
lscpu | grep "Model name" | sed 's/Model name: *//'
log_message "Total Memory:"
free -h | grep "Mem:" | awk '{print $2}'
log_message "Operating System:"
uname -a

# Load required modules
log_message "=== Loading Modules ==="
module load gcc/9.2.0
module load R/4.4.0
log_message "Loaded gcc/9.2.0 and R/4.4.0"

# Verify URD object exists
log_message "=== Checking Input Files ==="
if [ ! -f "data/initial_urd_object_20250210_1206.rds" ]; then
    log_message "Error: URD object file not found at data/initial_urd_object_20250210_1206.rds"
    exit 1
fi
log_message "Input URD object found"

# Log initial resource usage
log_message "=== Initial Resource Usage ==="
get_memory_usage
log_message "CPU Usage:"
top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1"%"}'

# Run dimensionality reduction analysis
log_message "=== Starting Analysis ==="
log_message "Running dimensionality reduction..."

# Time the R script execution
time Rscript run_dimensionality_reduction.R
SCRIPT_STATUS=$?

if [ $SCRIPT_STATUS -ne 0 ]; then
    log_message "Error in run_dimensionality_reduction.R (Exit code: $SCRIPT_STATUS)"
    exit 1
fi

# Calculate runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME-START_TIME))
HOURS=$((RUNTIME / 3600))
MINUTES=$(( (RUNTIME % 3600) / 60 ))
SECONDS=$(( (RUNTIME % 3600) % 60 ))

# Log final resource usage and completion
log_message "=== Final Resource Usage ==="
get_memory_usage
log_message "CPU Usage:"
top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1"%"}'

log_message "=== Analysis Summary ==="
log_message "Total Runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
log_message "Exit Status: $SCRIPT_STATUS"
log_message "Output files saved in results/"

# Check if expected output files exist
log_message "=== Output File Verification ==="
for file in "data/urd_object_with_dimred.rds" "data/urd_object_clean.rds" "results/parameter_summary.csv"; do
    check_file "$file"
done

log_message "Analysis completed"