#!/bin/bash
#SBATCH --job-name=SilvaDB    # Job name
#SBATCH --partition=normal              # Partition name (change to 'cpu' 360G or 'gpu' 740G as needed)
#SBATCH --mem=60G                    # Total memory per node (adjust as needed)
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=20                   # Number of tasks (usually 1 for single-job scripts)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --time=4-00:00:00            # Time limit (D-HH:MM:SS)
#SBATCH --output=/work/birdlab/databases/logs/silvaDB-%j.out     # Standard output and error log (%j will be replaced by job ID)

############################
### download_silvaDB.sh ###
############################

# This script downloads the most recent versions of the Silva Database

#----------------------------
# Execute:
# sbatch download_silvaDB.sh {locus}
#
# example:
# sbatch download_silvaDB.sh SSU
#----------------------------

locus=${1} #
#locus=SSU (or LSU)
database_dir=/work/birdlab/databases
script_dir=/work/birdlab/software/gcl_bioinformatic_tools/scripts

#0. Get URL for latest Silva ARB file
silva_url=https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_${locus}Ref_NR99_tax_silva.fasta.gz

#1 - Download the FASTA File
mkdir -p ${database_dir}/silva_latest/${locus}
cd ${database_dir}/silva_latest/${locus}
wget --no-check-certificate -c ${silva_url}
silva_dirty=$(realpath ${silva_url##*/})
silva_cleaned=${silva_dirty%.fasta.gz}_cleaned.fasta

#2. Download taxonomic mapping files
silva2ncbi_url=https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_${locus,,}_138.2.acc_taxid.gz
wget --no-check-certificate -c ${silva2ncbi_url}
silva2ncbi=$(realpath ${silva2ncbi_url##*/})

#3. FASTA: seqid -> accver
zcat "$silva_dirty" \
  | awk '
      /^>/{
        sub(/^>/,"",$1)
        seqid=$1
        n=split(seqid,p,"\\.")
        accver=p[1]"."p[2]      # accession.version
        print accver"\t"seqid
      }' \
  | sort -k1,1 \
  > accver2seqid.tsv

# 2) MAP: accver -> taxid
zcat "$silva2ncbi" \
  | awk 'NF>=2{
      id=$1; taxid=$2
      n=split(id,p,"\\.")
      accver=p[1]"."p[2]
      print accver"\t"taxid
    }' \
  | sort -k1,1 \
  > accver2taxid.tsv

# 3) Join: accver -> (seqid, taxid)  => seqid -> taxid  (BLAST taxid_map format)
join -t $'\t' -1 1 -2 1 accver2seqid.tsv accver2taxid.tsv \
  | awk -F'\t' '{print $2"\t"$3}' \
  > taxid_map

#4 - Clean the FASTA file
zcat ${silva_dirty} \
| awk '
  /^>/ {
    split($0,a," ")
    print a[1]
    next
  }
  { print }
' > ${silva_cleaned}

#5 - Format for RainbowBridge
BLAST=/work/birdlab/singularity_cache/ncbi-blast-latest.img
singularity exec --bind /work ${BLAST} \
  makeblastdb -in ${silva_cleaned} \
  -parse_seqids \
  -dbtype nucl \
  -taxid_map taxid_map \
  -out silva_latest

rm ${silva_cleaned} ${silva_dirty} accver2seqid.tsv accver2taxid.tsv ${silva2ncbi}
