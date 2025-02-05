import argparse
import os
from Bio import SeqIO

def extract_protein_from_genbank(genbank_files, protein_name, output_fasta):
    """
    Extracts a specific protein by name from multiple GenBank files and writes it to a single FASTA file.
    
    Parameters:
    genbank_files (list): List of paths to GenBank files.
    protein_name (str): Name of the protein to extract (e.g., "terminase large subunit").
    output_fasta (str): Path to the output FASTA file.
    """
    with open(output_fasta, "a") as fasta_out:  # Open in append mode
        for genbank_file in genbank_files:
            for record in SeqIO.parse(genbank_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS" and "product" in feature.qualifiers and "ID" in feature.qualifiers:
                        if protein_name in feature.qualifiers["product"]:
                            sequence = feature.qualifiers.get("translation", [""])[0]
                            if sequence:
                                protein_id = feature.qualifiers.get("ID", ["unknown"])[0]
                                header = f">{record.id}_{protein_id}_{protein_name}\n"
                                fasta_out.write(header + sequence + "\n")
    print(f"Extracted protein sequences saved to {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract marker proteins from GenBank files.")
    parser.add_argument("--genbank", nargs='+', required=True, help="Paths to GenBank files.")
    parser.add_argument("--protein", required=True, help="Protein name to extract.")
    parser.add_argument("--output", required=True, help="Output FASTA file.")
    
    args = parser.parse_args()
    extract_protein_from_genbank(args.genbank, args.protein, args.output)
