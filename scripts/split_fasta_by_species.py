#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Split a FASTA file by species/taxid. "
            "All sequences for a given species/taxid will go into one file."
        )
    )
    parser.add_argument("input_fasta", help="Path to the input FASTA file.")
    parser.add_argument(
        "output_dir",
        help="Directory where individual species-specific FASTA files will be written."
    )
    args = parser.parse_args()

    # 1) Read the FASTA and group records by species/taxid
    species_records = defaultdict(list)

    print("Reading and grouping sequences by species/taxID...")
    with open(args.input_fasta, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            # Example record.id: "MG559732.1.Clydonella_sawyeri_2201168"
            parts = record.id.split(".", 2)
            if len(parts) == 3:
                species_taxid = parts[2]  # e.g. "Clydonella_sawyeri_2201168"
            else:
                species_taxid = record.id  # fallback if format is unexpected

            species_records[species_taxid].append(record)

    # Create output directory if needed
    os.makedirs(args.output_dir, exist_ok=True)

    # 2) Write each group to its own FASTA file (one at a time)
    print("Writing grouped records to separate files...")
    for species_taxid, records in species_records.items():
        output_fasta_path = os.path.join(
            args.output_dir,
            f"{species_taxid}.fasta"
        )

        # Write all sequences for this species/taxid
        with open(output_fasta_path, "w") as out_fh:
            SeqIO.write(records, out_fh, "fasta")

        # Optional: print progress
        # print(f"Wrote {len(records)} sequences to {output_fasta_path}")

    print("Done. All files have been written.")

if __name__ == "__main__":
    main()

#Usage Example - made with chat GPT
#Use to split Cleaned MIDORI database by species to then identify outliers
#python split_fasta_by_species.py large_sequences.fasta output_species_dir
