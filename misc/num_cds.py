#!/usr/bin/env python
import argparse
import os
from collections import defaultdict
from Bio import SeqIO

def extract_genes(genbank_directory):
    genes_dict = defaultdict(lambda: defaultdict(int))
    file_names = []
    for file_name in os.listdir(genbank_directory):
        if file_name.endswith(".gbk") or file_name.endswith(".gb"):
            file_names.append(file_name)
            file_path = os.path.join(genbank_directory, file_name)
            for record in SeqIO.parse(file_path, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        gene_name = feature.qualifiers.get("product", ["Unknown_gene"])[0]
                        genes_dict[gene_name][file_name] += 1
    return genes_dict, file_names

def write_to_csv(genes_dict, file_names, output_file):
    with open(output_file, "w") as f:
        # Write header
        f.write("Gene," + ",".join(file_names) + "\n")
        # Write gene names and the number of CDS in each file
        for gene, counts_per_file in genes_dict.items():
            counts = [str(counts_per_file.get(file_name, 0)) for file_name in file_names]
            f.write(f"{gene},{','.join(counts)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script that takes a directory containing genbank files and writes gene presence absence table")
    parser.add_argument('-d', '--directory', dest='directory', help='Enter the directory containing the genbank files')
    parser.add_argument('-o', dest='output', help='Enter the output tabular format')
    results = parser.parse_args()
    genes_dict, file_names = extract_genes(results.directory)
    write_to_csv(genes_dict, file_names, results.output)
