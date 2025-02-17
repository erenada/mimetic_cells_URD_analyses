# Single-Cell Trajectory Analysis Pipeline

This repository contains R scripts and documentation for performing dimensionality reduction and trajectory analysis on single-cell RNA sequencing data using URD.

## Repository Structure

```mermaid
graph TD
    Root["/"] --> Doc1[README.md]
    Root --> Doc2[methodology_info.md]
    Root --> Scripts[scripts/]
    Root --> Test[test/]
    Root --> Data[data/]
    Root --> Results[results/]
    Root --> Logs[logs/]
    Root --> Resources[resources/]

    %% Scripts structure
    Scripts --> R[R/]
    Scripts --> Submit[submit/]
    
    R --> R1[01_create_urd_object.R]
    R --> R2[02_find_variable_genes.R]
    R --> R3[03_run_dimensionality_reduction.R]
    R --> R4[04_run_diffusion_map.R]
    R --> R5[05_run_pseudotime.R]

    Submit --> S1[01_submit_urd_object.sh]
    Submit --> S2[02_submit_var_genes.sh]
    Submit --> S3[03_submit_dimred.sh]
    Submit --> S4[04_submit_diffmap.sh]
    Submit --> S5[05_submit_pseudotime.sh]

    %% Test structure
    Test --> TestScripts[test_scripts/]
    Test --> TestData[test_data/]
    Test --> TestResults[test_results/]

    TestScripts --> T1[01_create_test_urd_object.R]
    TestScripts --> T2[02_run_test_variable_genes.R]
    TestScripts --> T3[03_run_test_dimensionality_reduction.R]
    TestScripts --> T4[04_run_test_diffusion_map.R]
    TestScripts --> T5[05_run_test_pseudotime.R]
    TestScripts --> T6[check_test_urd_object.R]

    %% Results structure
    Results --> Plots[plots/]
    Results --> VarGenes[variable_genes/]
    Results --> DimRed[dimensionality_reduction/]
    Results --> DiffMap[diffusion_map/]
    Results --> Pseudo[pseudotime/]

    Plots --> P1[dimensionality_reduction/]
    Plots --> P2[diffusion_map/]
    Plots --> P3[pseudotime/]
    Plots --> P4[variable_genes/]

    %% Styling
    classDef directory fill:#f9f,stroke:#333,stroke-width:2px
    classDef script fill:#bbf,stroke:#333,stroke-width:1px
    classDef doc fill:#fff,stroke:#333,stroke-width:1px
    
    class Scripts,R,Submit,Test,TestScripts,TestData,TestResults,Data,Results,Logs,Resources,Plots,VarGenes,DimRed,DiffMap,Pseudo,P1,P2,P3,P4 directory
    class R1,R2,R3,R4,R5,S1,S2,S3,S4,S5,T1,T2,T3,T4,T5,T6 script
    class Doc1,Doc2 doc
```

### Directory Overview

- 📁 **scripts/** - Main analysis scripts
  - 📁 **R/** - R analysis scripts
    - 📜 `01_create_urd_object.R` - Creates URD object
    - 📜 `02_find_variable_genes.R` - Identifies variable genes
    - 📜 `03_run_dimensionality_reduction.R` - PCA and tSNE
    - 📜 `04_run_diffusion_map.R` - Diffusion maps
    - 📜 `05_run_pseudotime.R` - Pseudotime calculation
  - 📁 **submit/** - O2 cluster submission scripts
    - 📜 `01_submit_urd_object.sh` → Submits URD creation
    - 📜 `02_submit_var_genes.sh` → Submits gene analysis
    - 📜 `03_submit_dimred.sh` → Submits dimension reduction
    - 📜 `04_submit_diffmap.sh` → Submits diffusion mapping
    - 📜 `05_submit_pseudotime.sh` → Submits pseudotime calc

- 📁 **test/** - Testing framework
  - 📁 **test_scripts/** - Test implementation
  - 📁 **test_data/** - Test datasets
  - 📁 **test_results/** - Test outputs

- 📁 **results/** - Analysis outputs
  - 📁 **plots/** - Visualizations
  - 📁 **variable_genes/** - Gene analysis results
  - 📁 **dimensionality_reduction/** - PCA/tSNE results
  - 📁 **diffusion_map/** - Diffusion mapping results
  - 📁 **pseudotime/** - Trajectory results

- 📁 **data/** - Input data (not tracked)
- 📁 **logs/** - Execution logs
- 📁 **resources/** - Additional documentation

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