#!/usr/bin/env python
from Bio import SeqIO
import os

def split_fasta(input_file, sample, output_dir):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Check if the input file is empty
    if os.path.getsize(input_file) == 0:
        # Touch the output file
        output_file = os.path.join(output_dir, f"{sample}_1.fasta")
        open(output_file, 'a').close()
        return

    # Read the input FASTA file
    with open(input_file, "r") as infile:
        count=1
        for record in SeqIO.parse(infile, "fasta"):
            # Construct the output file name
            output_file = os.path.join(output_dir, f"{sample}_{count}.fasta")

             # Rename the contig ID
            record.id = f"{sample}_{count}"
            record.description = ""  # Clear the description to avoid redundancy

            # Write the single record to the output file
            with open(output_file, "w") as outfile:
                SeqIO.write(record, outfile, "fasta")
            
            count  += 1

split_fasta(snakemake.input.fasta, snakemake.params.sample, snakemake.params.outdir)

# Example usage
#input_file = "path_to_your_input.fasta"
#output_dir = "path_to_output_directory"
#split_fasta(input_file, output_dir)
