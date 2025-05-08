# How to use DNA Quantification Shiny Apps

1. Login to gawain shiny server: [http://10.5.146.65/DNA_Quantification/](http://10.5.146.65/DNA_Quantification/) (must be in TAMUCC or using VPN)
2. Applications are designed to be run in order
	1. [`http://10.5.146.65/DNA_Quantification/1-quant_plate/`](http://10.5.146.65/DNA_Quantification/1-quant_plate/)
		- [`code/1 - quantPlate_shiny_2.R`](<code/1 - quantPlate_shiny_2.R>)
	2. [`http://10.5.146.65/DNA_Quantification/2-DNA_concentration/`](http://10.5.146.65/DNA_Quantification/2-DNA_concentration/)
	3. [`http://10.5.146.65/DNA_Quantification/3-DNA_normalization_PCR/`](http://10.5.146.65/DNA_Quantification/3-DNA_normalization_PCR)

# How to update Apps on Gawain

1. Copy updated shiny app to `gawain:/srv/shiny-server/DNA_Quantification/*/app.R`
	- Replace `*` with the appropriate app being updated
2. Logfiles are found here: `gawain:/var/log/shiny-server/`


# How to host apps on a different server

1. Setup Shiny-Server on new server: https://posit.co/products/open-source/shiny-server/
2. Ensure all used R packages are installed by the shiny-server user
3. Copy each shiny app into their own folders within the shiny-server