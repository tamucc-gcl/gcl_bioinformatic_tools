# How to use DNA Quantification Shiny Apps

1. Login to gawain shiny server: http://10.5.146.65/DNA_Quantification/ (must be in TAMUCC or using VPN)
2. 

# How to update Apps on Gawain

1. Copy updated shiny app to `gawain:/srv/shiny-server/DNA_Quantification/*/app.R`
	- Replace `*` with the appropriate app being updated
2. Logfiles are found here: `gawain:/var/log/shiny-server/`


# How to host apps on a different server

1. Setup Shiny-Server on new server: https://posit.co/products/open-source/shiny-server/
2. Ensure all used R packages are installed by the shiny-server user
3. Copy each shiny app into their own folders within the shiny-server