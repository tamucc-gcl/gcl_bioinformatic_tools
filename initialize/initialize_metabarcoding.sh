#Use to initialize the filetree and necessary scripts for generic metabarcoding project

master_dir=/scratch/group/p.bio240270.000/gcl_bioinformatic_tools
dirname=${1}


#Scripts
mkdir ${dirname}/scripts
cp ${master_dir}/scripts/run_rainbowBridge.sh ${dirname}/scripts
cp ${master_dir}/scripts/summarise_rainbowbridge.R ${dirname}/scripts

#Default YAML - need to edit as needed
cp ${master_dir}/initial_yamls/metabarcode_rainbowbridge_paired.yml ${dirname}
echo "Rename and fill in file-paths for YAML found ${dirname}"

#Build Folder Structure
mkdir ${dirname}/data
mkdir ${dirname}/output
cp ${master_dir}/initial_readmes/metabarcoding_readme.md ${dirname}/REAME.md

# Instructions for next steps
echo "Create metadata file if none exists"
echo "Populate Data"
echo "Personalize ReadME"
echo "Create as github repository"
echo "Create as DVC repository"
