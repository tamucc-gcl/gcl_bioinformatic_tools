git clone git@github.com:tamucc-gcl/gcl_bioinformatic_tools.git

# Make index
sudo cp gcl_bioinformatic_tools/labTools/DNA_Quant/code/index.html /srv/shiny-server/DNA_Quantification/

#Make app folders
sudo mkdir -p /srv/shiny-server/DNA_Quantification/1-quant_plate
sudo mkdir -p /srv/shiny-server/DNA_Quantification/2-DNA_concentration
sudo mkdir -p /srv/shiny-server/DNA_Quantification/3-DNA_normalization_PCR

#Copy App code
sudo cp gcl_bioinformatic_tools/labTools/DNA_Quant/code/1\ -\ quantPlate_shiny_v4.R  /srv/shiny-server/DNA_Quantification/1-quant_plate/app.R
sudo cp gcl_bioinformatic_tools/labTools/DNA_Quant/code/2\ -\ dna_amount_shiny_v10.R  /srv/shiny-server/DNA_Quantification/2-DNA_concentration/app.R
sudo cp gcl_bioinformatic_tools/labTools/DNA_Quant/code/3\ -\ Normalize_DNA_for_PCR_v4.R  /srv/shiny-server/DNA_Quantification/3-DNA_normalization_PCR/app.R

# Writeable folders for outputs
sudo mkdir -p /srv/shiny-server/DNA_Quantification/1-quant_plate/outdir \
              /srv/shiny-server/DNA_Quantification/2-DNA_concentration/outdir \
              /srv/shiny-server/DNA_Quantification/3-DNA_normalization_PCR/outdir
sudo chown -R shiny:shiny /srv/shiny-server/DNA_Quantification/1-quant_plate/outdir \
              /srv/shiny-server/DNA_Quantification/2-DNA_concentration/outdir \
              /srv/shiny-server/DNA_Quantification/3-DNA_normalization_PCR/outdir
sudo chmod -R u+rwx /srv/shiny-server/DNA_Quantification/1-quant_plate/outdir \
              /srv/shiny-server/DNA_Quantification/2-DNA_concentration/outdir \
              /srv/shiny-server/DNA_Quantification/3-DNA_normalization_PCR/outdir

# Extras for app 2
sudo cp gcl_bioinformatic_tools/labTools/DNA_Quant/code/dna_concentration*stan /srv/shiny-server/DNA_Quantification/2-DNA_concentration/model/
sudo mkdir -p /srv/shiny-server/DNA_Quantification/2-DNA_concentration/model
sudo mkdir -p /srv/shiny-server/DNA_Quantification/2-DNA_concentration/cmdstan_output
sudo chown -R shiny:shiny /srv/shiny-server/DNA_Quantification/2-DNA_concentration/cmdstan_output \
                          /srv/shiny-server/DNA_Quantification/2-DNA_concentration/model
sudo chmod -R u+rwx /srv/shiny-server/DNA_Quantification/2-DNA_concentration/cmdstan_output \
                    /srv/shiny-server/DNA_Quantification/2-DNA_concentration/model
#Check app 2 to ensure file locations are correct for model/output
