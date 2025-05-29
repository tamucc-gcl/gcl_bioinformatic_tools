# How to use DNA Quantification Shiny Apps

## 1  Login to the Shiny server

* URL <http://10.5.146.65/DNA_Quantification/>  
  *(on-campus or via VPN)*

## 2  Run the three apps **in order**

---

### 2.1  Quant-Plate

* **URL** <http://10.5.146.65/DNA_Quantification/1-quant_plate/>
* **Code** `code/1 - quantPlate_shiny_2.R`

#### Inputs

| Item | Example file / note |
|------|--------------------|
| Raw data file | `example_data/quant_rawData_accublue-nextgen.csv` |
| Plate-map file | `example_data/quant_plateMap_accublue-nextgen.csv` |
| Quant kit | *Autodetected* — “Accublue-nextgen” (pg) **or** “accublue / accuclear” (ng) |
| X variable | autofill → **rfu** |
| Y variable | autofill → **[pn]g_per_well** |
| Standard rows | autofill: rows with `plate_id == "standard"` |

#### Output

* Cleaned RFU table

#### Process (summary)

1. Detect kit → choose appropriate units  
2. Fill standards & variables automatically  
3. Export table for downstream modelling

---

### 2.2  DNA Concentration Modeller

* **URL** <http://10.5.146.65/DNA_Quantification/2-DNA_concentration/>
* **Code** `code/2 - dna_amount_shiny7.R`

#### Inputs

| Category | Details |
|----------|---------|
| **Data input** | Plate-map XLSX (`overall_dna-extract-plate-map.xlsx`) → choose sheet<br>Quant-plate output files (select **all**)<br>Stage: *original* or *requant* |
| **DNA concentrations** | Select column (autofill tries `input$y_var`) |
| **Model settings** | Chains · Iterations · Warm-up · Thinning |
| **Priors** (optional) | Set custom priors |
| **Flags** | Min pipettable volume (`input$min_volume`)<br>Max low-volume (`input$max_low_volume`)<br>Target DNA amount (`input$target_dna`)<br>“Excess DNA” multiplier (`input$mean_multiple`) |

#### Output

* CSV with model-based concentration estimates & dilution plan

#### Process (summary)

1. Merge plate-map and quant outputs  
2. Fit partial-pooling model to replicate data  
3. Calculate dilution factors, flag outliers & contaminants  
4. Produce per-sample recipe for normalisation

---

### 2.3  DNA Normalisation for PCR

* **URL** <http://10.5.146.65/DNA_Quantification/3-DNA_normalization_PCR/>
* **Code** `code/3 - Normalize_DNA_for_PCR.R`

#### Inputs

| Item | Description |
|------|-------------|
| Sample concentration CSV(s) | Output from App 2 |
| µL per PCR (`input$ul_per_PCR`) | Volume of DNA in PCR recipe |
| Number of PCR reactions (`input$number_PCR_rxns`) | Replicates you plan to run |
| DNA per PCR (`input$DNA_per_PCR`) | Desired ng of DNA in each reaction |
| Max volume per sample (`input$max_vol`) | Upper transfer-volume limit |

#### Output

* Transfer plate layout + water volumes  
* Per-reaction DNA mass table

#### Process (summary)

1. Compute total DNA needed (`rxns × DNA_per_PCR`)  
2. Derive goal volume given concentration  
3. Cap at `max_vol`; if exceeded, revise plan  
4. Output transfer + water volumes

---

## 3  Maintaining the Apps on **Gawain**

1. Copy updated `app.R` to  
   `gawain:/srv/shiny-server/DNA_Quantification/<app>/app.R`
2. Check logs at `/var/log/shiny-server/`

## 4  Hosting on another server

1. Install **Shiny-Server** – <https://posit.co/products/open-source/shiny-server/>
2. Install required R packages for the `shiny` user
3. Copy each app into its own folder inside Shiny-Server’s site dir
4. Copy this `index.html` (landing page) to  
   `/srv/shiny-server/DNA_Quantification/index.html`
