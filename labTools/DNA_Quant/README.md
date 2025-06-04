# DNA Quantification Toolkit

A comprehensive web-based toolkit for DNA quantification analysis, concentration calculation, and PCR normalization. Hosted on the TAMUCC Shiny server with three sequential applications designed for streamlined laboratory workflows.

## Quick Start

1. **Access the server:** [http://10.5.146.65/DNA_Quantification/](http://10.5.146.65/DNA_Quantification/) 
   > ⚠️ **Network Access Required:** Must be on TAMUCC network or connected via VPN

2. **Run applications in sequence:** Follow the numbered workflow below

3. **Use example data:** Sample files are provided for testing each application

## Workflow Overview

The toolkit consists of three sequential applications that process DNA quantification data from raw plate reader output to normalized PCR-ready volumes:

```
Raw Plate Data → [App 1] → Individual Concentrations → [App 2] → Mean Concentrations → [App 3] → PCR Volumes
```

---

## Application 1: Quantification Analysis

<details>
  <summary>Click to expand</summary>
   
**URL:** [`http://10.5.146.65/DNA_Quantification/1-quant_plate/`](http://10.5.146.65/DNA_Quantification/1-quant_plate/)

### Purpose
Process raw fluorescence data from plate readers and calculate DNA concentrations using optimized standard curves.

### Input Files
- **Raw Data File:** Fluorescence readings from plate reader
  - Example: [`quant_rawData_accublue-nextgen.csv`](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/example_data/quant_rawData_accublue-nextgen.csv)
- **Plate Map File:** Sample layout and standard information  
  - Example: [`quant_plateMap_accublue-nextgen.csv`](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/example_data/quant_plateMap_accublue-nextgen.csv)

### Configuration Parameters
| Parameter | Auto-Fill Behavior | Description |
|-----------|-------------------|-------------|
| **Quant Kit Used** | From filename | "Accublue-nextgen" (pg), "accublue", or "accuclear" (ng) |
| **X Variable** | "rfu" | Independent variable (fluorescence values) |
| **Y Variable** | "[pn]g_per_well" | Dependent variable (concentration units) |
| **Standard Rows** | From plate map | Wells with `plate_id` = "standard" |

### Process
1. Import fluorescence data and plate layout
2. Identify standard curve samples automatically
3. Fit optimal regression model to standards
4. Calculate unknown sample concentrations
5. Export individual replicate results

### Output
Individual DNA concentrations for each well/replicate, ready for statistical analysis in App 2.

</details>

---

## Application 2: Concentration Statistics

**URL:** [`http://10.5.146.65/DNA_Quantification/2-DNA_concentration/`](http://10.5.146.65/DNA_Quantification/2-DNA_concentration/)

### Purpose
Calculate mean DNA concentrations across replicates using Bayesian statistical models and prepare samples for PCR normalization.

### Input Files
- **Plate Map File:** Overall sample information
  - Example: [`overall_dna-extract-plate-map.xlsx`](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/example_data/overall_dna-extract-plate-map.xlsx)
- **Quant Output Files:** Results from Application 1 (select all relevant files)

### Key Parameters

#### DNA Analysis Settings
- **Excel Sheet Name:** Specify correct sheet in plate map file
- **DNA Concentration Column:** Auto-populated from App 1 output
- **Quant Stage:** Original quantification or re-quantification

#### Statistical Model Settings (Advanced)
- **Number of Chains:** Independent MCMC chains (default: 4)
- **Sampling Iterations:** Total posterior samples (default: 2000)
- **Warmup Iterations:** Burn-in samples (default: 1000)
- **Thinning Interval:** Sample spacing for independence (default: 1)

#### PCR Preparation Settings
| Parameter | Default | Description |
|-----------|---------|-------------|
| **Minimum Pipettable Volume** | 1 µL | Smallest volume for accurate pipetting |
| **Maximum Low Volume** | 5 µL | Max volume from low-concentration samples |
| **Target DNA Amount** | 10 ng/pg | Goal DNA quantity per reaction |
| **Excess DNA Threshold** | 3× | Flag samples above this multiple of mean |

### Statistical Process

#### Sample Grouping
The application automatically categorizes samples based on the `sample_type` column:
- **Samples:** Standard experimental samples
- **Controls:** Quality control samples (extraction, field, filter controls)

*Accepted `sample_type` values:* `'sample'`, `'control'`, `'extraction control'`, `'field control'`, `'filter control'`

#### Bayesian Model
Uses hierarchical modeling with partial pooling for robust concentration estimates:

```r
concentration ~ is_control + (1 | sample_id)
variance ~ is_control + (1 | sample_id)
```

#### Dilution Calculations
1. **Target template volume:** `(target_DNA × mean_concentration) / sample_concentration`
2. **Volume constraints:**
   - If target > max_volume → use max_volume
   - If target < min_volume → calculate dilution factor
3. **Dilution factor:** `sample_concentration / mean_concentration` (rounded for pipetting)
4. **Post-dilution metrics:** Adjusted concentrations and volumes

#### Quality Flags
- **Contaminated Controls:** Controls with DNA levels exceeding sample averages
- **Uncertain Estimates:** Samples with high concentration variance
- **Excess DNA:** Samples exceeding user-defined threshold above mean
- **High Variability:** Samples flagged for both uncertainty and excess DNA

---

## Application 3: PCR Normalization

**URL:** [`http://10.5.146.65/DNA_Quantification/3-DNA_normalization_PCR/`](http://10.5.146.65/DNA_Quantification/3-DNA_normalization_PCR/)

### Purpose
Calculate precise transfer volumes and water additions for normalized PCR reactions.

### Input Parameters
| Parameter | Description |
|-----------|-------------|
| **Sample Concentration CSV(s)** | Output files from Application 2 |
| **µL per PCR** | DNA volume used in each PCR reaction |
| **Number of PCR reactions** | Total reactions planned per sample |
| **DNA per PCR (ng)** | Target DNA amount per reaction |
| **Max volume per sample (µL)** | Maximum transfer volume limit |

### Calculations

#### Goal Volume Determination
```r
goal_volume = (number_PCR_reactions × DNA_per_PCR) / sample_concentration
```

#### Volume Optimization
```r
transfer_volume = min(goal_volume, max_volume)
water_to_add = (number_PCR_reactions × µL_per_PCR) - transfer_volume
actual_DNA_per_PCR = (sample_concentration × transfer_volume) / number_PCR_reactions
```

### Output
- **Transfer volumes:** Exact DNA volumes to pipette
- **Water additions:** Volume needed for normalization
- **Actual DNA amounts:** Real DNA quantity per reaction after normalization

---

## File Format Requirements

### Raw Data Files (App 1)
- **Format:** CSV with fluorescence values
- **Structure:** Wells in rows or columns with RFU measurements
- **Headers:** Include well identifiers and fluorescence readings

### Plate Map Files
- **Format:** CSV or Excel (.xlsx)
- **Required columns:**
  - Well identifiers matching raw data
  - `plate_id` (mark standards as "standard")
  - `sample_id` for sample tracking
  - `sample_type` for controls vs. samples

### Example File Structure
```
well_id,plate_id,sample_id,sample_type,rfu
A01,standard,std_1,standard,15000
A02,standard,std_2,standard,12000
B01,sample,fish_001,sample,8500
```

## Troubleshooting

### Common Issues

#### Application 1: Poor Standard Curves
- **Symptoms:** Low R² values, poor concentration estimates
- **Solutions:** 
  - Verify standard concentrations in plate map
  - Check for pipetting errors in standards
  - Remove outlier standard points if justified

#### Application 2: High Coefficient of Variation
- **Symptoms:** Large confidence intervals, uncertain flags
- **Solutions:**
  - Review pipetting technique for replicates
  - Check for systematic errors (edge effects, bubbles)
  - Consider re-quantifying problematic samples

#### Application 3: Unrealistic Volumes
- **Symptoms:** Very high or low transfer volumes
- **Solutions:**
  - Adjust target DNA amounts for your PCR requirements
  - Modify maximum volume constraints
  - Consider diluting high-concentration samples

### System Requirements
- **Browser:** Modern web browser with JavaScript enabled
- **Network:** TAMUCC network access or VPN connection
- **Files:** CSV and Excel file support

## Server Administration

### Updating Applications
1. Copy updated Shiny app to server location:
   ```bash
   scp app.R gawain:/srv/shiny-server/DNA_Quantification/[1-3]-*/app.R
   ```
2. Restart the Shiny server:
   ```bash
   sudo systemctl restart shiny-server
   ```

### Log Files
Server logs are located at: `/var/log/shiny-server/`

### Hosting on Different Servers
1. **Install Shiny Server:** Follow [Posit installation guide](https://posit.co/products/open-source/shiny-server/)
2. **Install R packages** for the `shiny-server` user
3. **Copy applications** to individual folders in `/srv/shiny-server/`
4. **Deploy landing page:** Copy `index.html` to `/srv/shiny-server/DNA_Quantification/`

## Source Code

- **Application 1:** [`code/1 - quantPlate_shiny_2.R`](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/code/1%20-%20quantPlate_shiny_2.R)
- **Application 2:** [`code/2 - dna_amount_shiny7.R`](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/code/2%20-%20dna_amount_shiny7.R)  
- **Application 3:** [`code/3 - Normalize_DNA_for_PCR.R`](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/code/3%20-%20Normalize_DNA_for_PCR.R)

## Support

For technical support or questions:
- **Primary Contact:** jason.selwyn@tamucc.edu
- **Institution:** Texas A&M University-Corpus Christi Genomics Core Lab
- **Documentation:** This README and inline application help

---

**Version:** 2.0  
**Last Updated:** 2025  
**Maintainer:** TAMUCC Genomics Core Lab
