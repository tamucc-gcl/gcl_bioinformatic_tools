#prj-name
PROJECT DESCRIPTION  - ADJUST

path to repo in Launch supercomputer `XXX`  - ADJUST

---

eDNA data has been processed with [rainbow_bridge](https://github.com/mhoban/rainbow_bridge) a fork of [eDNAFlow](https://github.com/mahsa-mousavi/eDNAFlow)

Repo Structure:

* data (Original Fastq sequence files with metadata and sample-mapping files - stored using [DVC](https://dvc.org/))
* [scripts](scripts) (Scripts and code used to analyze data)
* intermediate_files (modified fastq if any, along with rainbow_bridge working outputs - stored using [DVC](https://dvc.org/))
* [output](output) (Results)

All commands run in directory: [commands_run.md](commands_run.md)

---

## Data

DESCRIBE_DATA - ADJUST

Data files:

* [`data/*R[12]_[0-9]{3}.fastq.gz`]() - ADJUST

Parameter files:

* [metabarcode_rainbowbridge_paired.yml](metabarcode_rainbowbridge_paired.yml)
* [paired_demuxed_barcode.tsv](paired_demuxed_barcode.tsv)
* [sample.map](sample.map)
* [metadata.csv](metadata.csv)

## Protocol 
NOT DONE YET - TO EDIT

1. Files were downloaded from the TAMUCC grid via [`scripts/gridDownloader.sh`](scripts/gridDownloader.sh) using the link: [``]() - ADJUST

2. The integrity of the fastq format was then checked with [`scripts/checkFQ.sh`](scripts/checkFQ.sh) and [`scripts/validateFQPE.sbatch`](scripts/validateFQPE.sbatch)

3. Classified into zOTUs, classified to the lowest taxonomic level, and counted occurances in samples with [`scripts/run_rainbowBridge.sh`](scripts/run_rainbowBridge.sh)

	- Classifications are derived from the [MIDORI2 Unique Species database](https://www.reference-midori.info/index.html) and based on the lowest common ancestor (LCA) algorithm using the  which works by:
	
		```
		For each unique zOTU & each unique BLAST Hit (sequence ID) collapse to a single unique zOTU & BLAST Hit:
			- Selecting the first value for non-numeric fields (taxid to domain, length to slen, gapopen to stitle).
			- Taking the maximum value for pident, bitscore, qcov, and qcovhsp.
			- Taking the minimum value for mismatch and evalue.
		For each zOTU:
			- Count number of BLAST hits
			- For each Species:
				- Filter to the minimum e-value, followed by maximum pident, followed by maximum bitscore, followed by mininum mismatch
		Remove BLAST hits to zOTUs below the minimum query coverage and percent identity
			- Will remove zOTUs with no BLAST Hits above the thresholds
		For each zOTU: 		
			- Calculate percent identity difference between each blast hit and the highest percent identity BLAST hit within that zOTU
			- Remove hits with the difference greater than a maximum allowable difference
			- Find lowest shared taxonomic level across all remaining BLAST hits
		```
		
---

## Classification and Database Caveats 
* Incomplete Reference Sequences: The Midori2 database, like other reference databases, relies on curated sequences submitted by researchers. Some taxa, especially rare, cryptic, or poorly studied groups, may be underrepresented or entirely absent.

* Sequence Misannotations: Errors in the annotation of sequences (e.g., incorrect species or genus assignments) can propagate through analyses and lead to incorrect taxonomic identifications.

* Bias Toward Model or Economically Important Species: The database may be skewed toward species of greater economic or scientific interest, leaving gaps in representation for other taxa.

## Results
