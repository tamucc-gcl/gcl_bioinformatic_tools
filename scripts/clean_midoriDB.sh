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


#MASH=/scratch/group/p.bio240270.000/software/mash/mash
script_dir=/scratch/group/p.bio240270.000/gcl_bioinformatic_tools/scripts

in_fasta=${1}
min_sequences=${2}

# Check if stdout is a terminal; if yes, add the --bar option.
if [ -t 1 ]; then
    PARALLEL_OPTS="--bar -0 -j ${SLURM_CPUS_PER_TASK}"
else
    PARALLEL_OPTS="-0 -j ${SLURM_CPUS_PER_TASK}"
fi


#Split input fasta into species specific fasta files
module load Anaconda3; source activate biopython
python ${script_dir}/split_fasta_by_species.py ${in_fasta} output_species_dir
conda deactivate

#Align sequences
module load Anaconda3; source activate mafft
find output_species_dir -type f -name '*.fasta' -print0 \
| parallel ${PARALLEL_OPTS} --env min_sequences '
    count=$(grep -c "^>" "{}")
    if [ "$count" -lt "$min_sequences" ]; then
        exit 0
    fi
    mafft --thread 1 --quiet --auto "{}" > "{}.aligned"
'
conda deactivate

#Calculate pairwise distances between all sequences within all taxa
module load Anaconda3; source activate emboss
find output_species_dir -type f -name '*.aligned' -print0 \
| parallel ${PARALLEL_OPTS} --env min_sequences '
    count=$(grep -c "^>" "{}")
    if [ "$count" -lt "$min_sequences" ]; then
        exit 0
    fi
    distmat -sequence "{}" -nucmethod 1 -outfile "{}.lwrTri" 2>/dev/null
'
conda deactivate

#Identify groupings#
Rscript ${script_dir}/speciesDist_DirichletProcess.R output_species_dir
