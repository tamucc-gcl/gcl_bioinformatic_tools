# GCL Bioinformatics Tools

This repository contains general scripts, templates, and configurations used for GCL Bioinformatic projects.

## Folder Structure
- **initialize/**: Contains setup scripts for creating project directories and initializing workflows.
- **readmes/**: Templates and examples for project documentation.
- **scripts/**: Scripts for running and summarizing metabarcoding workflows.
- **yamls/**: Configuration files for the Rainbow Bridge pipeline.

## Getting Started
1. Use the initialization script to set up the project directory:
   ```bash
   bash initialize_metabarcoding.sh <project_directory>
   ```
2. Edit the YAML file copied from the `yamls` to the `<project_directory>` folder to reflect your data and parameters.
3. Populate new folder with data
4. Set-up git & dvc versioning
5. Perform analysis