#!/bin/bash
#SBATCH --job-name=DfamDB    # Job name
#SBATCH --partition=normal              # Partition name (change to 'cpu' 360G or 'gpu' 740G as needed)
#SBATCH --mem=60G                    # Total memory per node (adjust as needed)
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=20                   # Number of tasks (usually 1 for single-job scripts)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --time=4-00:00:00            # Time limit (D-HH:MM:SS)
#SBATCH --output=/work/birdlab/databases/logs/dfamDB-%j.out     # Standard output and error log (%j will be replaced by job ID)


############################
### download_midoriDB.sh ###
############################

# This script downloads the most recent versions of the Midori Curated Database

#----------------------------
# Execute:
# sbatch download_dfamDB.sh
#----------------------------

database_dir=/work/birdlab/databases
script_dir=/work/birdlab/software/gcl_bioinformatic_tools/scripts

mkdir -p ${database_dir}/dfam && cd dfam
wget -r -np -nH --cut-dirs=4 -A "dfam39_full.*.h5.gz,dfam39_full.*.h5.gz.md5,README.txt" \
  https://www.dfam.org/releases/current/families/FamDB/