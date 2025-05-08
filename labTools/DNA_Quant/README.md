# How to use DNA Quantification Shiny Apps

1. Login to gawain shiny server: [http://10.5.146.65/DNA_Quantification/](http://10.5.146.65/DNA_Quantification/) (must be in TAMUCC or using VPN)
2. Applications are designed to be run in order
	1. [`http://10.5.146.65/DNA_Quantification/1-quant_plate/`](http://10.5.146.65/DNA_Quantification/1-quant_plate/)
		- [`code/1 - quantPlate_shiny_2.R`](<code/1 - quantPlate_shiny_2.R>)
		- Inputs: 
		- Output:
		- Process: 
		
	2. [`http://10.5.146.65/DNA_Quantification/2-DNA_concentration/`](http://10.5.146.65/DNA_Quantification/2-DNA_concentration/)
		- [`code/2 - dna_amount_shiny7.R`](<code/2 - dna_amount_shiny7.R>)
		- Inputs: 
			-Data Input:
				- Plate Map File: 
				- Excel Sheet Name: Specify the sheet within the plate map file which is formatted properly (see example)
				- Upload Quant Output Files (select all): Output from App1
				- Quant Stage: Are these requants or original.
			- DNA Concentrations
				- Select DNA Concentration: Column containing DNA concentrations (attempts to auto-populate)
					- `input$y_var`
				- Model Settings: Optional model fitting controls (Advanced user setting)
					- Number of Chains: How many independent model chains to run
					- Number of Sampling Iterations: Total number of samples to draw from the posterior
					- Number of Warmup Iterations: How many of the total samples are used to parameterize fitting
					- Thinning Interval: How many samples to discard to get uncorrelated posterior samples
				- Priors Settings: Optional prior setting controls (Advanced user setting)
				- Flag Settings:
					- Minimum Pipettable Volume: The smallest volume that you would like to pipette (µL)
						- `input$min_volume`
					- Maximum Low Volume: The maximum volume to take from low concentration DNA samples (µL)
						- `input$max_low_volume `
					- Target Amount of DNA: The goal quantity of DNA to have in the normalized sample (ng or pg depending on input)
						- `input$target_dna`
					- Excess DNA is 'X' times more than the mean: Flag samples which are "X" times greater than the mean DNA concentration
						- `input$mean_multiple`
		- Output:
		- Process: 
			- Merges all the quantification output files with the plate map file
			- Based on quant replication samples calculates the average concentration of each sample
				- If there is a `sample_type` column splits estimates between samples and control (specified by including the word "control" in the sample_type column
					- Acceptable values in `sample_type` are: 'sample', 'control', 'extraction control', 'field control', and 'filter control'
				- Partial pooling of sample ID to improve individual estimates by using information from all samples
				- Variation in concentration allowed to vary by sample (ie. no assumption of equal variance)
				- `input$y_var ~ is_control + (1 | sample_id)` & `shape ~ is_control + (1 | sample_id)`
			- Calculate a dilution factor and post-dilution measures (if needed)
				- Calculate the target template volume to pipette
					- `ul_per_rxn = (input$target_dna * mean_dna_concentration) / !!sym(str_c(input$y_var, "_mean"))`
					- `mean_dna_concentration` is the calculated mean DNA concentration of all samples (excluding controls)
				- Volume limits:
					- If the target template volume is more than the the specified max volume then use the specified max volume
						- `ul_per_rxn > input$max_low_volume ~ input$max_low_volume`
					- If the target template volume is less than the minimum pipettable volume then calculate a dilution to perform
						- `ul_per_rxn < input$min_volume ~ NA_real_`
						- dilute to the overall average concentration of samples (excluding any controls)
							- `is.na(ul_per_rxn) ~ measured_conc / mean_dna_concentration`
							- The dilution factor is rounded for easier pipetting
						- Calculate the post-dilution concentration and volumes of DNA
							- `postDilution_concentration = if_else(dilution_factor > 0, measured_conc / dilution_factor, NA_real_)`
							- `postDilution_ul_per_rxn = round((input$target_dna * mean_dna_concentration) / postDilution_concentration)`
			- Calculate the amount of DNA for the reaction (post-dilution if indicated)
				- `rxn_ng = if_else(is.na(ul_per_rxn), postDilution_concentration * postDilution_ul_per_rxn, measured_conc * ul_per_rxn)`
			- Flag potential samples of interest
				- Any control sample where the upper 95% confidence interval of the DNA concentration is larger than the average of the non-control sample concentrations is flagged as a potential contaminant.
				- Any sample which has more dispersion than the 95% variance confidence interval is flagged as an uncertain estimate.
				- Any sample with a DNA concentration in excess of X times the mean DNA concentration (user specified)
				- Any sample which is both highly variable and as an excess of DNA
				
			
	3. [`http://10.5.146.65/DNA_Quantification/3-DNA_normalization_PCR/`](http://10.5.146.65/DNA_Quantification/3-DNA_normalization_PCR)
		- [`code/3 - Normalize_DNA_for_PCR.R`](<code/3 - Normalize_DNA_for_PCR.R>)
		- Inputs:
			- Upload Sample Concentration CSV(s)
				- Output from App2
			- µL per PCR: Volume of DNA used in the PCR recipe
				- `input$ul_per_PCR`
			- Number of PCR reactions: The number of PCR reactions you would like to be able to perform
				- `input$number_PCR_rxns`
			- DNA per PCR (ng): The amount of DNA desired in each PCR reaction
				- `input$DNA_per_PCR`
			- Max volume per sample (µL): The maximum volume of DNA you are willing to move to the transfer plate
				- `input$max_vol`
		- Outputs:
		- Process:
			- For each sample a goal volume of DNA is calculated based on the desired number of PCR reactions and the target quantity of DNA per pcr (user input) along with the DNA concentration of the sample.
				- `goal_volume_ul = (input$number_PCR_rxns * input$DNA_per_PCR) / ng_per_ul_mean`
			- If the goal volume of DNA required to do the specifed number of PCR reactions with the desired amount of DNA is greater than the specied maximum volume then default to using the maximum volume
				- `transfer_volume = if_else(goal_volume_ul > input$max_vol, input$max_vol, goal_volume_ul)`
			- Calculates the amount of water to be added to the transfer plate so each well can have the necessary DNA concentration to get the specified amount of DNA in each PCR reaction using the same volume of DNA.
				- `ul_to_add = (input$number_PCR_rxns * input$ul_per_PCR) - transfer_volume`
			- Reports the actual quantity of DNA in each PCR reaction:
				- `actual_ng_dna_per_pcr = (ng_per_ul_mean * transfer_volume) / input$number_PCR_rxns`

# How to update Apps on Gawain

1. Copy updated shiny app to `gawain:/srv/shiny-server/DNA_Quantification/*/app.R`
	- Replace `*` with the appropriate app being updated
2. Logfiles are found here: `gawain:/var/log/shiny-server/`


# How to host apps on a different server

1. Setup Shiny-Server on new server: https://posit.co/products/open-source/shiny-server/
2. Ensure all used R packages are installed by the shiny-server user
3. Copy each shiny app into their own folders within the shiny-server