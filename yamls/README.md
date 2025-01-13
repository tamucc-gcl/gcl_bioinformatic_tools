# YAMLs

This folder contains configuration files for the Rainbow Bridge workflow.

## Files
- **metabarcode_rainbowbridge_paired.yml**
  - YAML configuration file for running the Rainbow Bridge pipeline with paired-end metabarcoding data.
  - Key parameters:
    - Metadata paths
    - Barcode and sample-map paths
    - Blast database configuration
    - Workflow settings (e.g., quality thresholds, denoising, taxonomy collapsing).

## Usage
Modify the file paths and parameters in the YAML file as needed before running the workflow:
```yaml
metadata: PATH_TO_METADATA
barcode: PATH_TO_BARCODE
sample-map: PATH_TO_SAMPLEIDS
reads: PATH_TO_FASTQ
```