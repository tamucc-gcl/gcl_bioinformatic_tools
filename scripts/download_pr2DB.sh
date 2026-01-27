#!/bin/bash
#SBATCH --job-name=prrDB    # Job name
#SBATCH --partition=normal              # Partition name (change to 'cpu' 360G or 'gpu' 740G as needed)
#SBATCH --mem=60G                    # Total memory per node (adjust as needed)
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=20                   # Number of tasks (usually 1 for single-job scripts)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --time=4-00:00:00            # Time limit (D-HH:MM:SS)
#SBATCH --output=/work/birdlab/databases/logs/prrDB-%j.out     # Standard output and error log (%j will be replaced by job ID)
###############################################################################
# 0) HARDEN: consistent sorting & joins
###############################################################################
export LC_ALL=C

###############################################################################
# 1) User inputs (export before running, or edit here)
###############################################################################
LOCUS=SSU                              # SSU (18S), LSU, etc. (as PR2 provides)
PR2_VERSION=5.1.1               # PR2 release tag/version
DB_ROOT=/work/birdlab/databases/PR2_latest       # where you store databases
THREADS=${SLURM_CPUS_PER_TASK}

# PR2 FASTA (UTAX)
PR2_FASTA_URL=https://github.com/pr2database/pr2database/releases/download/v${PR2_VERSION}/pr2_version_${PR2_VERSION}_${LOCUS}_taxo_long.fasta.gz

# PR2 metadata Excel/TSV/CSV URL (YOU must set to the correct asset for the release)
# Examples (placeholder): https://.../pr2_version_5.1.1_metadata.xlsx
PR2_META_URL=https://github.com/pr2database/pr2database/releases/download/v${PR2_VERSION}/pr2_version_${PR2_VERSION}_merged.xlsx

# NCBI accession2taxid URLs
NCBI_GB_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
NCBI_WGS_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz"
NCBI_DEAD_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz"

###############################################################################
# 2) Output paths
###############################################################################
OUTDIR="${DB_ROOT}/${LOCUS}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

PR2_FASTA_GZ="$(basename "${PR2_FASTA_URL}")"
PR2_FASTA_CLEANED="${PR2_FASTA_GZ%.fasta.gz}_cleaned.fasta"

PR2_META_FILE="$(basename "${PR2_META_URL}")"
PR2_TAXO_TSV="pr2_accession__to__pr2_taxonomy.tsv"

PR2_TO_GB_TSV="pr2_accession__to__genbank_accession.tsv"  # pr2_seqid<TAB>key (accession-only / wgs stem / pseudo)
ACC2TAXID_TSV="key__to__taxid.tsv"                        # key<TAB>taxid
TAXID_MAP="taxid_map.tsv"                                 # pr2_seqid<TAB>taxid

# filtered outputs
PR2_FASTA_MAPPED="${PR2_FASTA_CLEANED%.fasta}_MAPPED_ONLY.fasta"
TAXID_MAP_MAPPED="${TAXID_MAP%.tsv}_MAPPED_ONLY.tsv"
PR2_TAXO_UNMAPPED="pr2_taxonomy_UNMAPPED_ONLY.tsv"

###############################################################################
# 3) Download PR2 FASTA
###############################################################################
echo "[1/9] Download PR2 FASTA"
wget -c -O "${PR2_FASTA_GZ}" "${PR2_FASTA_URL}"

###############################################################################
# 2. Download excel associated with that database version
###############################################################################
echo "[2/9] Download PR2 metadata: ${PR2_META_URL}"
wget -c -O "${PR2_META_FILE}" "${PR2_META_URL}"

###############################################################################
# 4) Clean FASTA headers (keep only seqid before first '|')
###############################################################################
echo "[3/9] Clean FASTA headers -> ${PR2_FASTA_CLEANED}"
zcat "${PR2_FASTA_GZ}" \
| awk -F'|' '
  /^>/{
    sub(/^>/,"",$1)
    print ">"$1
    next
  }
  {print}
' > "${PR2_FASTA_CLEANED}"

echo "  FASTA headers: $(grep -c '^>' "${PR2_FASTA_CLEANED}")"

###############################################################################
# 5) Extract PR2 taxonomy from taxo_long headers (seqid<TAB>rank1;rank2;...)
###############################################################################
echo "[4/9] Extract taxonomy -> ${PR2_TAXO_TSV}"
zcat "${PR2_FASTA_GZ}" \
| awk -F'|' '
  BEGIN{OFS="\t"}
  /^>/{
    sub(/^>/,"",$1)
    seqid=$1
    tax=""
    for(i=5;i<=NF;i++) tax = (tax=="" ? $i : tax ";" $i)
    print seqid, tax
  }
' > "${PR2_TAXO_TSV}"

###############################################################################
# 6) Parse PR2 metadata: PR2 accession + "GenBank accession" key
#    IMPORTANT: your observed metadata key is accession-only (AB012059), not AB012059.1
###############################################################################
echo "[5/9] Parse metadata Excel -> ${PR2_TO_GB_TSV}"

module load miniconda3
source activate r_env_pr2

Rscript --vanilla - "$PR2_META_FILE" "$PR2_TO_GB_TSV" <<'RSCRIPT'
args <- commandArgs(trailingOnly=TRUE)
meta_file <- args[1]
out_tsv   <- args[2]
message("meta_file = ", meta_file)
message("out_tsv   = ", out_tsv)

readxl::read_excel(meta_file, 
                   col_types = 'text') |>
  dplyr::distinct(pr2_accession, 
                  genbank_accession) |>
  readr::write_tsv(out_tsv, 
                   col_names = FALSE)
RSCRIPT
conda deactivate
echo "  PR2->key rows: $(wc -l < "${PR2_TO_GB_TSV}")"

###############################################################################
# 7) Build key->taxid using NCBI accession2taxid
#    - GB accessions (AB012059) map via nucl_gb (column 1)
#    - WGS stems (AAAA02002154) map via nucl_wgs (column 1)
#    - dead_nucl used to recover missing GB accessions (column 1)
###############################################################################
echo "[6/9] Download NCBI accession2taxid (if needed)"
wget -c -O nucl_gb.accession2taxid.gz   "${NCBI_GB_URL}"
wget -c -O nucl_wgs.accession2taxid.gz  "${NCBI_WGS_URL}"
wget -c -O dead_nucl.accession2taxid.gz "${NCBI_DEAD_URL}"

cut -f2 "${PR2_TO_GB_TSV}" | awk 'NF{print $1}' | sort -u > pr2_keys.all
echo "  unique metadata keys: $(wc -l < pr2_keys.all)"

# Split keys
grep -E '^[A-Z]{1,2}[0-9]{5,8}$' pr2_keys.all > pr2_keys.gb || true
grep -E '^[A-Z]{4}[0-9]{8}$'     pr2_keys.all > pr2_keys.wgs || true
grep -Ev '(^[A-Z]{1,2}[0-9]{5,8}$|^[A-Z]{4}[0-9]{8}$)' pr2_keys.all > pr2_keys.other || true

echo "  keys: GB=$(wc -l < pr2_keys.gb) WGS=$(wc -l < pr2_keys.wgs) OTHER=$(wc -l < pr2_keys.other)"

extract_col1() {
  local gz="$1"
  local keep="$2"
  zcat "$gz" | awk -F'\t' 'NR==FNR{keep[$1]=1; next} NR>1 && keep[$1]{print $1"\t"$3}' "$keep" -
}

echo "  mapping GB keys from nucl_gb (col1)..."
extract_col1 nucl_gb.accession2taxid.gz pr2_keys.gb | awk '!seen[$1]++' > acc2taxid.gb.tsv

# dead_nucl recovery for GB
cut -f1 acc2taxid.gb.tsv | sort -u > pr2_keys.gb.mapped || true
comm -23 pr2_keys.gb pr2_keys.gb.mapped > pr2_keys.gb.missing || true
echo "  GB keys missing after nucl_gb: $(wc -l < pr2_keys.gb.missing)"

echo "  mapping remaining GB keys from dead_nucl (col1)..."
extract_col1 dead_nucl.accession2taxid.gz pr2_keys.gb.missing | awk '!seen[$1]++' > acc2taxid.dead.tsv

echo "  mapping WGS stems from nucl_wgs (col1)..."
extract_col1 nucl_wgs.accession2taxid.gz pr2_keys.wgs | awk '!seen[$1]++' > acc2taxid.wgs.tsv

cat acc2taxid.gb.tsv acc2taxid.dead.tsv acc2taxid.wgs.tsv | awk '!seen[$1]++' > "${ACC2TAXID_TSV}"
echo "  mapped key->taxid: $(wc -l < "${ACC2TAXID_TSV}") / $(wc -l < pr2_keys.all)"

###############################################################################
# 8) Build TAXID_MAP (PR2 accession -> taxid) by joining PR2_TO_GB_TSV with ACC2TAXID_TSV
###############################################################################
echo "[7/9] Build taxid_map -> ${TAXID_MAP}"

sort -k2,2 "${PR2_TO_GB_TSV}" > pr2_to_key.sorted.tsv
sort -k1,1 "${ACC2TAXID_TSV}" > key2taxid.sorted.tsv

# join output columns are: key  pr2_seqid  taxid
join -t $'\t' -1 2 -2 1 pr2_to_key.sorted.tsv key2taxid.sorted.tsv \
| awk -F'\t' '{print $2"\t"$3}' \
| awk '!seen[$1]++' \
> "${TAXID_MAP}"

echo "  taxid_map rows: $(wc -l < "${TAXID_MAP}")"
echo "  TAXID_MAP first 5:"
head -n 5 "${TAXID_MAP}" | sed -n 'l'

echo "FASTA first 3:"
grep '^>' "${PR2_FASTA_CLEANED}" | head -n 3 | sed -n 'l'

echo "TAXID_MAP first 3:"
head -n 3 "${TAXID_MAP}" | sed -n 'l'

# should now show overlap
join -t $'\t' -1 1 -2 1 \
  <(grep '^>' "${PR2_FASTA_CLEANED}" | sed 's/^>//' | awk '{print $1}' | sort -u | head -n 2000) \
  <(cut -f1 "${TAXID_MAP}" | sort -u | head -n 2000) \
| head

# Checks
comm -23 <(cut -f1 "${TAXID_MAP}" | sort -u) \
         <(grep '^>' "${PR2_FASTA_CLEANED}" | sed 's/^>//' | awk '{print $1}' | sort -u) \
| head

total=$(grep -c '^>' "${PR2_FASTA_CLEANED}")
mapped=$(wc -l < "${TAXID_MAP}")
echo "total=$total mapped=$mapped missing=$((total-mapped))"


###############################################################################
# 9) Build KEEP/UNMAPPED lists safely (NO comm pitfalls)
#    We will use join for set ops.
###############################################################################
echo "[8/9] Build mapped-only FASTA and mapped-only TAXID_MAP"

# All seqids from FASTA
grep '^>' "${PR2_FASTA_CLEANED}" | sed 's/^>//' | awk '{print $1}' | sort -u > pr2_seqids.all
# Mapped seqids from TAXID_MAP
cut -f1 "${TAXID_MAP}" | sort -u > pr2_seqids.mapped

# Compute unmapped = all - mapped using join (robust)
# join outputs intersection; so unmapped is those NOT in mapped:
# Use awk keep array instead (fast + robust), but we keep locale C anyway.
awk 'NR==FNR{m[$1]=1; next} !($1 in m){print $1}' pr2_seqids.mapped pr2_seqids.all > pr2_seqids.unmapped
echo "  seqids: total=$(wc -l < pr2_seqids.all) mapped=$(wc -l < pr2_seqids.mapped) unmapped=$(wc -l < pr2_seqids.unmapped)"

# Sanity check
if [[ $(( $(wc -l < pr2_seqids.mapped) + $(wc -l < pr2_seqids.unmapped) )) -ne $(wc -l < pr2_seqids.all) ]]; then
  echo "ERROR: mapped + unmapped != total (set logic failed). Aborting."
  exit 1
fi

# Filter FASTA to mapped seqids
awk '
  BEGIN{ while((getline < "pr2_seqids.mapped")>0) keep[$1]=1; close("pr2_seqids.mapped") }
  /^>/{ id=$0; sub(/^>/,"",id); split(id,a," "); seq=a[1]; p=(seq in keep) }
  { if(p) print }
' "${PR2_FASTA_CLEANED}" > "${PR2_FASTA_MAPPED}"

# Filter TAXID_MAP to mapped seqids (paranoia-safe)
join -t $'\t' -1 1 -2 1 \
  <(sort -k1,1 "${TAXID_MAP}") \
  <(sort -k1,1 pr2_seqids.mapped) \
| awk -F'\t' '{print $1"\t"$2}' \
> "${TAXID_MAP_MAPPED}"

echo "  mapped-only FASTA headers: $(grep -c '^>' "${PR2_FASTA_MAPPED}")"
echo "  mapped-only taxid_map rows: $(wc -l < "${TAXID_MAP_MAPPED}")"

if [[ $(grep -c '^>' "${PR2_FASTA_MAPPED}") -ne $(wc -l < "${TAXID_MAP_MAPPED}") ]]; then
  echo "ERROR: FASTA headers != TAXID_MAP rows after filtering. Aborting."
  exit 1
fi

# Create taxonomy TSV for unmapped seqids so you can inspect/drop them
sort -k1,1 "${PR2_TAXO_TSV}" > pr2_taxo.sorted.tsv
sort -u pr2_seqids.unmapped > pr2_seqids.unmapped.sorted
join -t $'\t' -1 1 -2 1 pr2_seqids.unmapped.sorted pr2_taxo.sorted.tsv > "${PR2_TAXO_UNMAPPED}"

echo "  wrote unmapped taxonomy TSV: ${PR2_TAXO_UNMAPPED} (rows: $(wc -l < "${PR2_TAXO_UNMAPPED}"))"

echo "[9/9] DONE"
echo "  Use these downstream:"
echo "    FASTA:     ${PR2_FASTA_MAPPED}"
echo "    TAXID_MAP: ${TAXID_MAP_MAPPED}"
echo "  Inspect unmapped taxonomy:"
echo "    ${PR2_TAXO_UNMAPPED}"

# Make BLAST DB
# Requires BLAST+ installed and available in PATH (makeblastdb)
BLAST=/work/birdlab/singularity_cache/ncbi-blast-latest.img
mv ${TAXID_MAP_MAPPED} taxid_map
singularity exec --bind /work ${BLAST} \
  makeblastdb \
    -in "${PR2_FASTA_MAPPED}" \
    -dbtype nucl \
    -parse_seqids \
    -taxid_map taxid_map \
    -out pr2_latest

rm ${PR2_FASTA_GZ} ${PR2_FASTA_CLEANED} ${PR2_META_FILE} ${PR2_TO_GB_TSV} ${ACC2TAXID_TSV} ${TAXID_MAP}
rm nucl_gb.accession2taxid.gz nucl_wgs.accession2taxid.gz dead_nucl.accession2taxid.gz