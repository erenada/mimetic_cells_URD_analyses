# Single-Cell Trajectory Analysis Pipeline

This repository contains R scripts and documentation for performing dimensionality reduction and trajectory analysis on single-cell RNA sequencing data using URD.

## Repository Structure

```
.
├── README.md                    # Main documentation
├── methodology_info.md          # Detailed methodology documentation
├── scripts/                     # Main analysis scripts
│   ├── R/                      # R analysis scripts
│   │   ├── 01_create_urd_object.R
│   │   ├── 02_find_variable_genes.R
│   │   ├── 03_run_dimensionality_reduction.R
│   │   ├── 04_run_diffusion_map.R
│   │   └── 05_run_pseudotime.R
│   └── submit/                 # O2 submission scripts
│       ├── 01_submit_urd_object.sh
│       ├── 02_submit_var_genes.sh
│       ├── 03_submit_dimred.sh
│       ├── 04_submit_diffmap.sh
│       └── 05_submit_pseudotime.sh
├── test/                       # Test-related files and scripts
│   ├── test_scripts/          # Test R scripts
│   │   ├── 01_create_test_urd_object.R
│   │   ├── 02_run_test_variable_genes.R
│   │   ├── 03_run_test_dimensionality_reduction.R
│   │   ├── 04_run_test_diffusion_map.R
│   │   ├── 05_run_test_pseudotime.R
│   │   └── check_test_urd_object.R
│   ├── test_data/             # Test dataset files
│   └── test_results/          # Test analysis results
├── data/                       # Main data directory (not tracked)
├── results/                    # Main analysis results
│   ├── plots/
│   │   ├── dimensionality_reduction/
│   │   ├── diffusion_map/
│   │   ├── pseudotime/
│   │   └── variable_genes/
│   ├── dimensionality_reduction/
│   ├── diffusion_map/
│   ├── pseudotime/
│   └── variable_genes/
├── logs/                       # Log files directory
└── resources/                  # Additional resources and documentation
```

## Scripts

### Analysis Scripts (in `scripts/R/`)
1. `01_create_urd_object.R`: Creates a URD object
2. `02_find_variable_genes.R`: Identifies variable genes across stages
3. `03_run_dimensionality_reduction.R`: Performs PCA, tSNE, and clustering analyses
4. `04_run_diffusion_map.R`: Calculates diffusion map components
5. `05_run_pseudotime.R`: Computes pseudotime ordering

### Submission Scripts (in `scripts/submit/`)
1. `01_submit_urd_object.sh`: Submits URD object creation
2. `02_submit_var_genes.sh`: Submits variable genes analysis
3. `03_submit_dimred.sh`: Submits dimensionality reduction
4. `04_submit_diffmap.sh`: Submits diffusion map calculation
5. `05_submit_pseudotime.sh`: Submits pseudotime analysis

## Documentation

- `methodology_info.md`: Comprehensive documentation of the methodology, including:
  - Parameter selection rationale
  - Data-driven parameter adjustments
  - Visualizations of parameter relationships
  - Literature references

## Requirements

- R >= 4.2.2
- URD package
- Seurat package
- RColorBrewer package

## Usage

### Local Execution
Run the scripts in the following order:
```R
source("scripts/R/01_create_urd_object.R")
source("scripts/R/02_find_variable_genes.R")
source("scripts/R/03_run_dimensionality_reduction.R")
source("scripts/R/04_run_diffusion_map.R")
source("scripts/R/05_run_pseudotime.R")
```

### O2 Cluster Execution
Submit jobs in the following order:
```bash
cd scripts/submit
sbatch 01_submit_urd_object.sh
# Wait for completion, then:
sbatch 02_submit_var_genes.sh
# Wait for completion, then:
sbatch 03_submit_dimred.sh
# Wait for completion, then:
sbatch 04_submit_diffmap.sh
# Wait for completion, then:
sbatch 05_submit_pseudotime.sh
```

## Results

The analysis generates:
1. URD object creation
2. Variable genes analysis
   - Variable genes list
   - Stage-specific statistics
3. Dimensionality reduction
   - PCA plots and statistics
   - tSNE visualizations
   - Clustering results
4. Diffusion map analysis
   - Diffusion components
   - Quality metrics
5. Pseudotime analysis
   - Pseudotime ordering
   - Stage progression plots
   - Stability assessment

Results are saved in the `results/` directory with the following structure:
```
results/
├── plots/
│   ├── dimensionality_reduction/
│   ├── diffusion_map/
│   └── pseudotime/
├── variable_genes/
├── dimensionality_reduction/
├── diffusion_map/
├── pseudotime/
└── variable_genes/
```

## References

See `methodology_info.md` for a comprehensive list of references and citations. 