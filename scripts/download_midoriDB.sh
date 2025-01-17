#!/bin/bash
#SBATCH --job-name=MidoriDB    # Job name
#SBATCH --partition=cpu              # Partition name (change to 'cpu' 360G or 'gpu' 740G as needed)
#SBATCH --mem=340G                    # Total memory per node (adjust as needed)
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=192                   # Number of tasks (usually 1 for single-job scripts)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --time=2-00:00:00            # Time limit (D-HH:MM:SS)
#SBATCH --output=midoriDB-%j.out     # Standard output and error log (%j will be replaced by job ID)


############################
### download_midoriDB.sh ###
############################

# This script downloads the most recent versions of the Midori Curated Database

#----------------------------
# Execute:
# sbatch download_midoriDB.sh
#
# example:
# sbatch download_midoriDB.sh
#----------------------------

script_dir=/scratch/group/p.bio240270.000/gcl_bioinformatic_tools/scripts

#1 - Identify most recently added fasta on MIDORI
module load Anaconda3; source activate r_env
Rscript ${script_dir}/get_midori_url.R
