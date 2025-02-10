# URD Package Manual

## Overview
URD (named after the Norse mythological figure) is an R package designed for reconstructing transcriptional trajectories from single-cell RNA-sequencing data. It specializes in modeling developmental processes as branching trees, making it particularly useful for studying cell specification and differentiation.

## Key Features
- Reconstruction of developmental trajectories
- Visualization of branching cellular decisions
- Analysis of gene expression changes along developmental paths
- Support for large-scale single-cell RNA-seq datasets
- Doublet detection using NMF module analysis
- Pseudotime analysis capabilities
- Advanced visualization tools for data inspection

## System Requirements

### Software Prerequisites
1. R (tested on versions 3.5 and 3.6)
2. RStudio
3. X11 window client (XQuartz for Mac users)
4. udunits2 library (pre-installed on Windows/Mac, requires manual installation on Linux)

### Hardware Requirements
- Sufficient RAM for large-scale single-cell data analysis
- Graphics capability for 3D visualization

## Installation

### Quick Installation
```R
source("https://raw.githubusercontent.com/farrellja/URD/master/URD-Install.R")
```

### Manual Installation Steps
1. Install devtools: `install.packages("devtools")`
2. Install Bioconductor dependencies:
```R
source("https://bioconductor.org/biocLite.R")
biocLite(c('Biobase', 'S4Vectors', 'AnnotationDbi', 'destiny'))
```
3. Install optional dependencies:
```R
biocLite(c('sva', 'rhdf5', 'scran'))
```
4. Install URD:
```R
library(devtools)
install_github(repo="farrellja/URD")
```

## Workflow Steps

### 1. Data Preparation
- Import single-cell RNA-seq data
- Quality control and filtering
- Normalization
- Feature selection

### 2. Trajectory Analysis
- Create URD object
- Calculate diffusion map
- Define root cells
- Build transition matrix
- Construct developmental trajectory
- Identify branch points

### 3. Visualization
- Plot trajectory tree
- Visualize gene expression
- Generate pseudotime plots
- Create branch point visualizations

### 4. Advanced Analysis
- Module detection
- Doublet removal
- Gene expression analysis relative to pseudotime
- Branch point analysis

## Available Tutorials

1. **Quick Start Tutorial**
   - Focus on axial mesoderm development
   - Available in both R Markdown and MD formats
   - Includes example dataset

2. **Supplementary Analyses**
   - Zebrafish embryogenesis (3.3-12 hours post-fertilization)
   - Adult Hydra development

## Example Datasets
1. Zebrafish Embryogenesis Dataset
   - Available from Broad Single-cell Portal
   - Covers 3.3-12 hours post-fertilization

2. Adult Hydra Dataset
   - Available from both Broad Single-cell Portal and Data Dryad
   - Demonstrates stem cell differentiation trajectories

## Version Information
- Current Version: 1.1.1 (as of April 28, 2020)
- Major Updates:
  - Version 1.1: Added NMF module analysis, doublet detection, and new visualization tools
  - Version 1.0: Initial release with core trajectory reconstruction functionality

## Support and Resources
- GitHub Repository: https://github.com/farrellja/URD
- Published Papers:
  1. "Single-cell reconstruction of developmental trajectories during zebrafish embryogenesis" (Science, 2018)
  2. "Stem cell differentiation trajectories in Hydra resolved at single-cell resolution" (Science, 2019)

## Troubleshooting Tips
1. DLL Error Resolution:
   - Increase R_MAX_NUM_DLLS to 200 in .Renviron
   - Command for Mac: `echo "R_MAX_NUM_DLLS=200" >> ~/.Renviron`

2. Installation Issues:
   - Consider using a separate conda environment
   - Check udunits2 installation
   - Verify all dependencies are properly installed

3. Permission Issues:
   - Specify custom library location when installing on clusters
   - Use appropriate user permissions

This manual provides a comprehensive overview of the URD package. For specific analyses or advanced usage, refer to the tutorials and example datasets provided in the repository. 