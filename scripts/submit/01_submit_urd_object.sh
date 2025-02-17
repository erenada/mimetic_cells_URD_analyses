#!/bin/bash

#SBATCH -J URD_create                       # Job name
#SBATCH -c 20                               # Request 20 cores
#SBATCH -t 6:00:00                         # Runtime in D-HH:MM format
#SBATCH -p short                            # Partition to run in
#SBATCH --mem=32G                           # Memory total in MiB (for all cores)
#SBATCH -o logs/create_urd_%j.out           # File to which STDOUT will be written, %j is job ID
#SBATCH -e logs/create_urd_%j.err           # File to which STDERR will be written, %j is job ID
#SBATCH --mail-type=ALL                     # Type of email notification- ALL=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=eren_ada@hms.harvard.edu
#SBATCH --job-name=URD_create              # Job name

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to get memory usage
get_memory_usage() {
    free -h | awk 'NR==2{printf "Memory Usage: %s/%s (%.2f%%)\n", $3,$2,$3*100/$2 }'
}

# Create necessary directories
mkdir -p logs data results

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

# Log initial resource usage
log_message "=== Initial Resource Usage ==="
get_memory_usage
log_message "CPU Usage:"
top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1"%"}'

# Run URD object creation
log_message "=== Starting Analysis ==="
log_message "Running URD object creation..."

# Time the R script execution
time Rscript scripts/R/01_create_urd_object.R
SCRIPT_STATUS=$?

if [ $SCRIPT_STATUS -ne 0 ]; then
    log_message "Error in create_urd_object.R (Exit code: $SCRIPT_STATUS)"
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
log_message "Output files saved in data/"

# Check if expected output files exist
log_message "=== Output File Verification ==="
if [ -f "data/initial_urd_object.rds" ]; then
    FILE_SIZE=$(stat --format=%s "data/initial_urd_object.rds")
    FILE_SIZE_MB=$(echo "scale=2; $FILE_SIZE/1024/1024" | bc)
    log_message "✓ data/initial_urd_object.rds successfully created ($FILE_SIZE_MB MB)"
else
    log_message "✗ data/initial_urd_object.rds not found"
    exit 1
fi

log_message "Analysis completed successfully" 