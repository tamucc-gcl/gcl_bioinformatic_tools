#!/bin/bash
#SBATCH --job-name=cleanMidori    # Job name
#SBATCH --partition=cpu              # Partition name (change to 'cpu' 360G or 'gpu' 740G as needed)
#SBATCH --mem=340G                    # Total memory per node (adjust as needed)
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=192                   # Number of tasks (usually 1 for single-job scripts)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --time=2-00:00:00            # Time limit (D-HH:MM:SS)
#SBATCH --output=midoriDB-%j.out     # Standard output and error log (%j will be replaced by job ID)


############################
### clean_midoriDB.sh ###
############################

# This script downloads the most recent versions of the Midori Curated Database

#----------------------------
# Execute:
# sbatch clean_midoriDB.sh <fasta>
#
# example:
# sbatch clean_midoriDB.sh MIDORI2_UNIQ_NUC_GB263_CO1_RAW_cleaned.fasta
#----------------------------


MASH=/scratch/group/p.bio240270.000/software/mash/mash
script_dir=/scratch/group/p.bio240270.000/gcl_bioinformatic_tools/scripts

in_fasta=${1}

#Split input fasta into species specific fasta files
module load Anaconda3; source activate biopython
python ${script_dir}/split_fasta_by_species.py ${in_fasta} output_species_dir
conda deactivate

#Align sequences


#Calculate pairwise distances between all sequences within all taxa
find output_species_dir -type f -name '*.fasta' -print0 \
| parallel -0 -j "${SLURM_CPUS_PER_TASK}" "${MASH} triangle -E -p 1 {} > {}.lwrTriangle"

Rscript ${script_dir}/speciesDist_DirichletProcess.R output_species_dir
