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
midori_url=$(Rscript ${script_dir}/get_midori_url.R)
echo "Downloading: ${midori_url}"
conda deactivate

#2 - Download the FASTA File
mkdir -p /scratch/group/p.bio240270.000/databases/midori2_latest
cd /scratch/group/p.bio240270.000/databases/midori2_latest
wget -c ${midori_url}
midori_dirty=$(realpath ${midori_url##*/})
midori_cleaned=${midori_dirty%.fasta.gz}_cleaned.fasta

#3 - Clean-up database
# Cleanning includes:
# (a) Keeping only the accession number, species name (or hybrid), and species' NCBI taxonomic ID (removes extra taxonomic info such as higher taxo-levels as well as lower level such as subpecies, stratins,etc)
	#However, if lower level info did exit for a record (subpecies, strains, etc). You rainbow_bridge will give taxonomic info at this level but higher levels can then be deduce from this by the user
# (b) Truncates names to 50 characters

# Cleaning fasta
zcat "${midori_dirty}" | \
	# Truncates first column to the first comma
	awk -F'\t' '{ split($1, a, ","); $1 = a[1]; print }' | \
	# adds a ; at the end of name lines for further processing
	sed '/^>/ s/$/;/'  | \
	# gets rids of white space and everything after species info
	sed -e 's/\(;species_[^;]*_[0-9]*\);.*/\1/'  -e 's/ /_/g' | \
	# keeps only accession number and species info
	awk '/^>/ {split($0, a, ";species_"); split(a[1], b, "."); print b[1] "." b[2] "." a[2]} !/^>/ {print}' | \
	# truncate names to 49 characters
	awk '/^>/ {print substr($0, 1, 47)} !/^>/ {print}' | \
	# add a serial number at the end '-1', '-2', etc to duplicated names.
	awk '/^>/ {count[$0]++; if (count[$0] > 1) $0 = $0 "-" count[$0]-1} 1' > "${midori_cleaned}"


#4 - Create taxid_map
# column1= seq names (makeblastdb has a 50 characters limit) and
# column2= ncbi taxonomic ID (this will be taken from the midori2 raw fasta in case there were issues while "cleaning the fasta  in the previous step"

grep '>' "${midori_cleaned}" | sed 's/>//' > column1
zgrep '>' "${midori_dirty}" | sed 's/.*_//' > column2
paste column1 column2 > taxid_map

# remove temp files
rm column1 column2 ${midori_dirty}

#Clean-up the database by removing sequences too different from their species ID
#bash ${script_dir}/clean_midoriDB.sh "${midori_cleaned}"

#5 - Format for RainbowBridge
BLAST=/scratch/group/p.bio240270.000/software/ncbi-blast-latest.img
singularity exec --bind /scratch,/ztank ${BLAST} \
  makeblastdb -in ${midori_cleaned} \
  -parse_seqids \
  -dbtype nucl \
  -taxid_map taxid_map \
  -out midori2_latest

rm ${midori_cleaned}
