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
				- Model Settings: Optional model fitting controls (Advanced user setting)
					- Number of Chains: How many independent model chains to run
					- Number of Sampling Iterations: Total number of samples to draw from the posterior
					- Number of Warmup Iterations: How many of the total samples are used to parameterize fitting
					- Thinning Interval: How many samples to discard to get uncorrelated posterior samples
				- Priors Settings: Optional prior setting controls (Advanced user setting)
					- Intercept Prior Mean: 
					- Intercept Prior SD: 
					- Beta Prior Mean: 
					- Beta Prior SD: 
					- Variance Prior Parameter 1: 
					- Variance Prior Parameter 2: 
					- Intercept Shape Prior Mean: 
					- Intercept Shape Prior SD: 
					- Beta Shape Prior Mean: 
					- Beta Shape Prior SD: 
					- Variance Shape Prior Parameter 1: 
					- Variance Shape Prior Parameter 2: 
				- Flag Settings:
					- Minimum Pipettable Volume: 
					- Maximum Low Volume:
					- Target Amount of DNA: 
		- Output:
		- Process: 
			- Merges all the quantification output files with the plate map file
			- Based on quant replication calculates 
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