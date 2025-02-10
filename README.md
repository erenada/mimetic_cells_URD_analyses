# Single-Cell Trajectory Analysis Pipeline

This repository contains R scripts and documentation for performing dimensionality reduction and trajectory analysis on single-cell RNA sequencing data using URD.

## Repository Structure

```
.
├── README.md
├── methodology_info.md          # Detailed methodology documentation
├── run_dimensionality_reduction.R  # Script for dimensionality reduction
├── run_linage_analysis.R       # Script for lineage analysis
├── create_urd_object.R         # Script for creating URD object
├── data/                       # Directory for data (not tracked)
└── results/                    # Directory for analysis results
    ├── plots/
    │   └── dimensionality_reduction/
    └── variable_genes/
```

## Scripts

1. `create_urd_object.R`: Creates URD object from Seurat object
2. `run_linage_analysis.R`: Performs lineage analysis and calculates variable genes
3. `run_dimensionality_reduction.R`: Performs PCA, tSNE, and clustering analyses

## Documentation

- `methodology_info.md`: Comprehensive documentation of the methodology, including:
  - Parameter selection rationale
  - Data-driven parameter adjustments
  - Visualizations of parameter relationships
  - Literature references

## Requirements

- R >= 4.0.0
- URD package
- Seurat package
- RColorBrewer package

## Usage

1. Place your input data in the `data/` directory
2. Run the scripts in the following order:
   ```R
   source("create_urd_object.R")
   source("run_linage_analysis.R")
   source("run_dimensionality_reduction.R")
   ```

## Results

The analysis generates:
1. Dimensionality reduction plots
2. Variable genes lists
3. Parameter selection summaries
4. Outlier detection results

Results are saved in the `results/` directory.

## References

See `methodology_info.md` for a comprehensive list of references and citations. 