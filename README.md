# Predicting_vector_abundance_in_Lithuania

Using hierarchical spatial-temporal models implemented in R-INLA to understand and predict patterns of insect vector abundance in Lithuanian lakes and rivers

R-code to reproduce analysis in Baker & Palinauskas (n.d.). Ecological integrity drives increased abundances of insect vectors (Culicidae, Simuliidae, Ceratopogonidae). Submitted to Parasites & Vectors

# Diptera Taxonomic Analysis Project

## Overview
This repository contains the complete analysis pipeline for analyzing Diptera taxonomic indices across different environmental conditions. The analysis consists of 8 sequential R scripts that must be run in order due to data dependencies.

## Prerequisites
- R (version ≥ 4.0.0)
- RStudio
- Required R packages (listed in each script)
- Access to external datasets (see Data Requirements below)

## Repository Structure
```
├── Data/
│   ├── 1_raw_macroinvertebrate_data_long.csv
│   ├── 2_ecological_ratios_and_classes.csv
│   └── 3_haemosporidian_parasite_data_long.csv
├── Outputs/
├── Plots/
├── Sensitivity/
├── Additional data/ (hidden in .gitignore)
├── Corine2018/ (user must download)
├── Scripts/
│   ├── 1_calculating_taxonomic_indices.R
│   ├── 2_extracting_corine_landcover_2018.R
│   ├── 3_extracting_terraclimate.R
│   ├── 4_extracting_elevation.R
│   ├── 5_model_implementation.R
│   ├── 6_model_predictions.R
│   ├── 7_plotting_parasite_data.R
│   └── 8_creating_covariate_panel_plot.R
└── README.md
```

## Data Requirements
Some datasets are too large for GitHub and must be downloaded separately:

1. **Corine Landcover 2018**
   - Download from: [URL to be provided]
   - Place in: `Corine2018/` directory
   
2. **TerraClimate Data**
   - Download instructions provided in script: `3_extracting_terraclimate.R`
   - Place in: `Additional data/TerraClimate/`

3. **Lithuanian Rivers Shapefile** (optional)
   - Request access from Lithuanian Environmental Protection Agency
   - Place in: `Additional data/GeoDatabase/`

## Analysis Workflow

### Step 1: Calculate Taxonomic Indices
**Script:** `1_calculating_taxonomic_indices.R`
- **Input:** 
  - `Data/1_raw_macroinvertebrate_data_long.csv`
  - `Data/2_ecological_ratios_and_classes.csv`
- **Output:** `Outputs/1_diptera_taxonomic_indices.csv`

### Step 2: Extract Corine Landcover Data
**Script:** `2_extracting_corine_landcover_2018.R`
- **Input:**
  - `Outputs/1_diptera_taxonomic_indices.csv`
  - `Corine2018/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif`
  - `Additional data/Corine Landcover/clc_legend.csv`
  - `Additional data/Corine Landcover/CLC2018_CLC2018_V2018_20_QGIS.txt`
- **Output:** `Outputs/2_diptera_taxonomic_indices_wCorine2018.csv`
- **Warning:** User must download Corine Landcover data

### Step 3: Extract TerraClimate Data
**Script:** `3_extracting_terraclimate.R`
- **Input:**
  - `Outputs/2_diptera_taxonomic_indices_wCorine2018.csv`
  - `Additional data/TerraClimate/linked_terraclimate_data.RDS`
- **Output:** `Outputs/3_diptera_taxonomic_indices_wCorine2018_TerraClimate.csv`
- **Warning:** User must download TerraClimate data

### Step 4: Extract Elevation Data
**Script:** `4_extracting_elevation.R`
- **Input:** `Outputs/3_diptera_taxonomic_indices_wCorine2018_TerraClimate.csv`
- **Output:** `Outputs/4_diptera_taxonomic_indices_wCorine2018_TerraClimate_elevation.csv`

### Step 5: Implement Model
**Script:** `5_model_implementation.R`
- **Input:** `Outputs/4_diptera_taxonomic_indices_wCorine2018_TerraClimate_elevation.csv`
- **Output:** 
  - `Outputs/5_unique_sites_for_plotting.csv`
  - Main plots: `Figure4`, `Figure6`
  - Supplement plots: `TableS2`, `FigureS1`, `FigureS6`, `FigureS7`
  - Sensitivity plots: Multiple figures (S2-S8)
- **Warning:** If `rgeoboundaries` fails, use `rnaturalearth` package

### Step 6: Generate Model Predictions
**Script:** `6_model_predictions.R`
- **Input:** `Outputs/4_diptera_taxonomic_indices_wCorine2018_TerraClimate_elevation.csv`
- **Output:** `Plots/Figure5_predicted_fixed_effects_without_spatial.png`

### Step 7: Plot Parasite Data
**Script:** `7_plotting_parasite_data.R`
- **Input:** `Data/3_haemosporidian_parasite_data_long.csv`
- **Output:** 
  - `Plots/Figure1_parasite_prevalence_dynamics.png`
  - `Plots/Figure1_parasite_prevalence_dynamics.RDS`

### Step 8: Create Covariate Panel Plot
**Script:** `8_creating_covariate_panel_plot.R`
- **Input:**
  - `Outputs/5_unique_sites_for_plotting.csv`
  - `Additional data/GeoDatabase/UETK_2024-05-02.gdb`
  - Corine and TerraClimate data
- **Output:** 
  - `Plots/Figure2_Sampling_sites_wWater.png`
  - `Plots/Figure3_covariate_panel_plot.png`
- **Warning:** Lithuanian river shapefile optional; TerraClimate required

## Workflow Diagram

```mermaid
graph TD
    subgraph Input Files
        R[Raw Data - macroinvertebrate_data.csv, ecological_ratios.csv]
        P[Parasite Data - haemosporidian_data.csv]
        E[External Data - Corine, TerraClimate, GeoDatabase]
    end
    
    subgraph Analysis Pipeline
        R --> S1[1. Taxonomic Indices - Creates taxonomic metrics]
        S1 --> S2[2. Corine Landcover - Adds landcover data]
        E --> S2
        S2 --> S3[3. TerraClimate - Adds climate variables]
        E --> S3
        S3 --> S4[4. Elevation - Adds elevation data]
        S4 --> S5[5. Model Implementation - Spatial regression models]
        S5 --> S6[6. Model Predictions - Generate predictions]
        S5 --> Sites[Unique Sites File]
        Sites --> S8[8. Covariate Panel - Create final visualizations]
        E --> S8
        P --> S7[7. Parasite Plots - Prevalence analysis]
    end
    
    subgraph Output Figures
        S5 --> O1[Model Results - Figure 4, 6, Table S2, Figures S1-S8]
        S6 --> O2[Predictions - Figure 5]
        S7 --> O3[Parasite Dynamics - Figure 1]
        S8 --> O4[Site Maps and Panels - Figure 2, 3]
    end
    
    %% Styling
    classDef input fill:#e1f5fe,stroke:#0288d1,stroke-width:2px;
    classDef script fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px;
    classDef output fill:#e8f5e9,stroke:#388e3c,stroke-width:2px;
    classDef intermediate fill:#fff3e0,stroke:#f57c00,stroke-width:2px;
    
    class R,P,E input;
    class S1,S2,S3,S4,S5,S6,S7,S8 script;
    class O1,O2,O3,O4 output;
    class Sites intermediate;
```

## Running the Analysis
1. Clone this repository
2. Download required external datasets (see Data Requirements)
3. Run scripts sequentially from 1 to 8
4. Output plots and files will be generated in respective directories

## Troubleshooting
- **Package Loading Issues:** If `rgeoboundaries` fails, use `rnaturalearth` as alternative
- **Missing Data:** Ensure all external datasets are downloaded and placed in correct directories
- **Memory Issues:** Some operations may require significant memory; close other applications if needed

## Citation
[Add citation information here]

## Contact
[Add contact information]

## License
[Add license information]
