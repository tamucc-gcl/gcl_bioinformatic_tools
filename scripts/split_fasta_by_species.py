#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description="Split a FASTA file by species/taxid on-the-fly."
    )
    parser.add_argument("input_fasta", help="Path to the input FASTA file.")
    parser.add_argument(
        "output_dir",
        help="Directory where individual species-specific FASTA files will be written."
    )
    args = parser.parse_args()

    # Create the output directory if it doesn't already exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Dictionary to hold open file handles keyed by species/taxid
    file_handles = {}

    # Parse the input FASTA and write each record immediately
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        # Example record.id = "MG559732.1.Clydonella_sawyeri_2201168"
        # Split on '.' only twice to handle sequences like "MG559732.1.some_name"
        parts = record.id.split('.', 2)
        if len(parts) == 3:
            # The third part after splitting on '.' is "Clydonella_sawyeri_2201168"
            species_taxid = parts[2]
        else:
            # Fallback if the format is unexpected
            species_taxid = record.id

        # Construct output FASTA path
        output_fasta_path = os.path.join(args.output_dir, f"{species_taxid}.fasta")

        # Open a new file handle if this species hasn't been encountered yet
        if species_taxid not in file_handles:
            file_handles[species_taxid] = open(output_fasta_path, "w")

        # Write the record to the appropriate file in FASTA format
        # You can either manually format or use SeqIO.write with a single-element list:
        SeqIO.write([record], file_handles[species_taxid], "fasta")

    # Close all file handles
    for fh in file_handles.values():
        fh.close()

if __name__ == "__main__":
    main()

#Usage Example - made with chat GPT
#Use to split Cleaned MIDORI database by species to then identify outliers
#python split_fasta_by_species.py large_sequences.fasta output_species_dir
