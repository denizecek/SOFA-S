# SOFA-S (Survival Omics Factor Analysis – Spatial) 
R-based pipeline for identifying spatial latent factors from single-cell spatial omics data using *MOFA2/MEFISTO* and evaluating their association with patient survival.

## Features

- Spatial binning of single-cell coordinates
- Aggregation of marker expression within spatial bins
- MOFA2/MEFISTO latent factor analysis
- Spatial visualization of latent factors
- Factor weight and variance explained plots
- Kaplan–Meier survival analysis
- Cox proportional hazards modelling
- ROI-level and bin-level factor summaries
- Cell-type integration and spatial overlay of dominant cell types
- Automated generation of publication-ready figures

## Input

The pipeline requires three input tables:

- **Coordinate table** – Cell coordinates (X, Y) and ROI information.
- **Omics table** – Single-cell marker expression matrix.
- **Metadata table** – Sample metadata including ROI, survival time and event status.

## Example

```r
results <- run_sofa(
  coord_table_path = "...",
  omic_data_path = "...",
  metadata_path = "...",
  output_dir = "...",
  bin_size = 500
)
```

## Notes
SOFA-S is under active development.
