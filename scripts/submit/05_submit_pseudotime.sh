#!/bin/bash

#SBATCH -J URD_pseudotime                   # Job name
#SBATCH -c 20                               # Request 20 cores
#SBATCH -t 12:00:00                         # Runtime in D-HH:MM format
#SBATCH -p short                            # Partition to run in
#SBATCH --mem=64G                           # Memory total in MiB (for all cores)
#SBATCH -o logs/pseudotime_%j.out           # File to which STDOUT will be written, %j is job ID
#SBATCH -e logs/pseudotime_%j.err           # File to which STDERR will be written, %j is job ID
#SBATCH --mail-type=ALL                     # Type of email notification- ALL=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=eren_ada@hms.harvard.edu

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to get memory usage
get_memory_usage() {
    free -h | awk 'NR==2{printf "Memory Usage: %s/%s (%.2f%%)\n", $3,$2,$3*100/$2 }'
}

# Create necessary directories
mkdir -p logs data results/pseudotime results/plots/pseudotime

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

# Verify input file exists
log_message "=== Checking Input Files ==="
if [ ! -f "data/urd_object_with_dm.rds" ]; then
    log_message "Error: URD object with diffusion map not found at data/urd_object_with_dm.rds"
    log_message "Please run 04_submit_diffmap.sh first"
    exit 1
fi
log_message "Input URD object found"

# Log initial resource usage
log_message "=== Initial Resource Usage ==="
get_memory_usage
log_message "CPU Usage:"
top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1"%"}'

# Run pseudotime analysis
log_message "=== Starting Analysis ==="
log_message "Running pseudotime calculation..."

# Time the R script execution
time Rscript scripts/R/05_run_pseudotime.R
SCRIPT_STATUS=$?

if [ $SCRIPT_STATUS -ne 0 ]; then
    log_message "Error in run_pseudotime.R (Exit code: $SCRIPT_STATUS)"
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
log_message "Output files saved in results/pseudotime/"

# Check if expected output files exist
log_message "=== Output File Verification ==="
if [ -f "data/urd_object_with_pseudotime.rds" ]; then
    FILE_SIZE=$(stat --format=%s "data/urd_object_with_pseudotime.rds")
    FILE_SIZE_MB=$(echo "scale=2; $FILE_SIZE/1024/1024" | bc)
    log_message "✓ data/urd_object_with_pseudotime.rds successfully created ($FILE_SIZE_MB MB)"
else
    log_message "✗ data/urd_object_with_pseudotime.rds not found"
    exit 1
fi

if [ -f "results/pseudotime/parameters.csv" ]; then
    log_message "✓ Parameters file created"
else
    log_message "✗ Parameters file not found"
fi

if [ -f "results/pseudotime/stage_statistics.csv" ]; then
    log_message "✓ Stage statistics file created"
else
    log_message "✗ Stage statistics file not found"
fi

# Check for plot files
for plot in "stage_pseudotime_distribution.png" "pseudotime_diffusion_components.png" "pseudotime_stability.png"; do
    if [ -f "results/plots/pseudotime/${plot}" ]; then
        log_message "✓ ${plot} created"
    else
        log_message "✗ ${plot} not found"
    fi
done

# Add periodic memory monitoring
monitor_resources() {
    while true; do
        log_message "=== Resource Monitor ==="
        get_memory_usage
        log_message "CPU Usage:"
        top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1"%"}'
        sleep 300  # Check every 5 minutes
    done
}

# Start resource monitoring in background
monitor_resources &
MONITOR_PID=$!

# After R script completion, stop the monitoring
kill $MONITOR_PID

log_message "Analysis completed successfully" 