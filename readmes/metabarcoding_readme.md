#prj-name
PROJECT DESCRIPTION  - ADJUST

path to repo in Launch supercomputer `XXX`  - ADJUST

---

eDNA data has been processed with [rainbow_bridge](https://github.com/mhoban/rainbow_bridge) a fork of [eDNAFlow](https://github.com/mahsa-mousavi/eDNAFlow)

Repo Structure:

* data (Original Fastq sequence files with metadata and sample-mapping files)
* scripts (Scripts and code used to analyze data)
* intermediate_files (modified fastq if any, along with rainbow_bridge working outputs)
* output (Results)

---

## Data

DESCRIBE_DATA - ADJUST

Data files:

* data/*R[12]_[0-9]{3}.fastq.gz
* [data/paired_demuxed_barcode.tsv](data/paired_demuxed_barcode.tsv)
* [data/sample.map](data/sample.map)
* [data/metadata.csv](data/metadata.csv)

Parameter files:

* [metabarcode_rainbowbridge_paired.yml](metabarcode_rainbowbridge_paired.yml)
