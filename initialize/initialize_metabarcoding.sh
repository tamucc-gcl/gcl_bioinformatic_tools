#Use to initialize the filetree and necessary scripts for generic metabarcoding project

master_dir=/scratch/group/p.bio240270.000/gcl_bioinformatic_tools
dirname=${1}


#Scripts
#mkdir -p ${dirname}/scripts
#cp ${master_dir}/scripts/run_rainbowBridge.sh ${dirname}/scripts
#cp ${master_dir}/scripts/summarise_rainbowbridge.R ${dirname}/scripts

#Default YAML - need to edit as needed
cp ${master_dir}/yamls/metabarcode_rainbowbridge_paired.yml ${dirname}
echo "Rename and fill in file-paths for YAML found ${dirname}"

#Build Folder Structure
mkdir -p ${dirname}/data
mkdir -p ${dirname}/intermediate_files
mkdir -p ${dirname}/output
cp ${master_dir}/readmes/metabarcoding_readme.md ${dirname}/README.md
touch ${dirname}/commands_run.md

# Instructions for next steps
echo "Any code run copy the lines to : ${dirname}/commands_run.md"
echo "Create metadata file if none exists"
echo "Populate Data"
echo "Update YAML to correct filepaths"
echo "Personalize ReadME - changes where readme says 'ADJUST'"
echo "Create as github repository"
echo "Create as DVC repository"

echo ""
echo ""
echo ""
echo "After above run by cd to ${dirname} and use:"
#echo "sbatch scripts/run_rainbowBridge.sh metabarcode_rainbowbridge_paired.yml intermediate_files"
echo "sbatch /work/birdlab/software/rainbow_bridge/run_rainbowBridge.sbatch metabarcode_rainbowbridge_paired.yml intermediate_files"
echo "Prior to running the first time you will need to build the rainbow bridge summarization conda environment"
echo "See this file for instructions /work/birdlab/software/rainbow_bridge/README.md"
echo ""
echo ""
echo ""
