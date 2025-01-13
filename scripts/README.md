# Scripts

This folder contains scripts used for running GCL Bioinformatic projects and summarizing results.

## Metabarcoding Files
- **run_rainbowBridge.sh**
  - SLURM job script to run the Rainbow Bridge workflow for metabarcoding analyses.
  - Usage: `sbatch run_rainbowBridge.sh <params YAML file> <output_dir>`

- **summarise_rainbowbridge.R**
  - R script for summarizing and visualizing outputs from the Rainbow Bridge workflow.
  - Includes functionality for creating summary metrics, taxonomic visualizations, and sample composition analyses.

## Notes
Ensure that the correct dependencies and modules (e.g., `GCC`, `rainbow_bridge`, `R`) are loaded in the computing environment before running these scripts.

