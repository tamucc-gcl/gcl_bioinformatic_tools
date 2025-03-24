#!/bin/bash
#SBATCH --job-name=rainbow_bridge    # Job name
#SBATCH --partition=cpu              # Partition name (change to 'cpu' 360G or 'gpu' 740G as needed)
#SBATCH --mem=340G                    # Total memory per node (adjust as needed)
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=192                   # Number of tasks (usually 1 for single-job scripts)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --time=2-00:00:00            # Time limit (D-HH:MM:SS)
#SBATCH --output=rainbow_bridge-%j.out     # Standard output and error log (%j will be replaced by job ID)
#SBATCH --mail-type=END,FAIL              #Send email on all job events
#SBATCH --mail-user=@tamucc.edu    #Send all emails to email_address
##SBATCH --ntasks=8                   #Request 8 tasks
##SBATCH --ntasks-per-node=2          #Request 2 tasks/cores per node

############################
### run_rainbow_bride.sh ###
############################

# This script runs rainbow_bridge in the TAMU Supercomputer 'Launch'. After logging in, you must enter a computing node in Launch.
# You will need a reference database, and a barcode decode file if your sequences have not been demultiplex

#----------------------------
# Execute:
# sbatch run_rainbow_brigde.sh <params YAML file> <output_dir> final
#
# example:
# sbatch run_rainbow_brigde.sh paired_unmuxed.yml output
# OR to run finalization copying use:
# sbatch run_rainbow_brigde.sh paired_unmuxed.yml output final
#----------------------------

# note: I am using NCBI's nucleotide database blastn v2_16_0. You might want to check if a newer version exist in https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# currently this database lives in '/home/u.eg195763/GCL/databases/ncbi_2_16_0/nt'
# "nt" is the name of the NCBI's nucleotide database (without file extensions)

# Load necessary modules (Launch - TAMUCC)
module load GCC/13.2.0 WebProxy rainbow_bridge/2024.07.15

# Print SLURM configuration
#env | grep SLURM

# Record the start time
start_time=$(date +%s)
echo "Start time: $(date)"
echo ""

# Get the name and full path of the script executed
#script_path=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $3}')
script_path=$(scontrol show job $SLURM_JOB_ID | awk -F'=' '/Command=/{print $2}')
script_dir=${script_path%/*}
echo ${script_dir} #TEST

# Initial reporting. Print parameters used
PARAMSFILE=$1
outdir=${2}
final=${3}


# Write the script name to the output file
echo "The script executed is: $script_path"
echo "Script executed from: $(pwd)"
echo -e "Using params file=$PARAMSFILE\n"
echo "cat $PARAMSFILE"
cat $PARAMSFILE
echo ""


mkdir -p ${outdir}
cp ${PARAMSFILE} ${outdir}/
cd ${outdir}

# Run rainbow_bridge with desired options
nextflow run -params-file $PARAMSFILE /scratch/group/p.bio240270.000/software/rainbow_bridge/rainbow_bridge.nf

#Summarize and produce various figures of RainbowBridge Outputs
module load Anaconda3; source activate r_env
Rscript ${script_dir}/summarise_rainbowbridge.R $(pwd) ${PARAMSFILE}
conda deactivate


## Copy files and finalize outputs to return to customer ##
if [ "${final}" == "final" ]; then
  echo "Finalizing Outputs"

  #Move outputs from intermediate to outputs
  mkdir -p ../output/rainbow_bridge

  mv *html ../output/rainbow_bridge/
  mv *png ../output/rainbow_bridge/
  mv *.csv ../output/rainbow_bridge/
  mv *metadata.txt ../output/rainbow_bridge/
  mv *fasta ../output/rainbow_bridge/
  cp -Lr output/phyloseq/phyloseq.rds ../output/rainbow_bridge/


  #If multiqc was not run because only one fastq file run multiqc for the sake of consistent output
  if ! test -f output/fastqc/initial/multiqc_report.html; then
    module load Anaconda3; source activate multiqc
    multiqc output/fastqc/initial/ -o output/fastqc/initial/
    conda deactivate
  fi

  if ! test -f output/fastqc/filtered/multiqc_report.html; then
    module load Anaconda3; source activate multiqc
    multiqc output/fastqc/filtered/ -o output/fastqc/filtered/
    conda deactivate
  fi

  mkdir -p ../output/fastqc/initial
  mkdir -p ../output/fastqc/filtered
  cp -Lr output/fastqc/initial/*html ../output/fastqc/initial/
  cp -Lr output/fastqc/filtered/*html ../output/fastqc/filtered/

  ## Tarball the full rainbow bridge analysis
  mkdir -p ../output/rainbowBridge_complete_output
  cp -Lr output/* ../output/rainbowBridge_complete_output
  cd ../output
  tar cvJf rainbow_bridge/rainbowBridge_complete_output.tar.xz rainbowBridge_complete_output
  rm -rf rainbowBridge_complete_output

  #Copy to file sharing server


fi


# Record the end time
end_time=$(date +%s)
echo ""
echo "End time: $(date)"

# Calculate the elapsed time
elapsed_time=$((end_time - start_time))

# Convert elapsed time to hours, minutes, and seconds
hours=$((elapsed_time / 3600))
minutes=$(( (elapsed_time % 3600) / 60 ))
seconds=$((elapsed_time % 60))

echo -e "Total time taken: ${hours}h ${minutes}m ${seconds}s\n"

echo -e "\nReport Resource Usage:"
seff $SLURM_JOB_ID
