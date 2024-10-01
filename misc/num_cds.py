#!/usr/bin/env python
import argparse
import os
from collections import defaultdict
from Bio import SeqIO
import re
import csv


#Code from https://github.com/linsalrob/EdwardsLab/blob/e49085a1b0c97735f93bc1a1261514b4829e0ef3/roblib/functions.py#L10-L59
def is_hypothetical(func):
    """
    Returns True if the function is hypothetical. Otherwise returns false
    :param func: string
    :return: boolean
    """

    if not func: return True
    if func.lower() == 'hypothetical protein': return True
    if re.search(r'lmo\d+ protein', func, re.IGNORECASE): return True
    if re.search(r'hypoth', func, re.IGNORECASE): return True
    if re.search(r'conserved protein', func, re.IGNORECASE): return True
    if re.search(r'gene product', func, re.IGNORECASE): return True
    if re.search(r'interpro', func, re.IGNORECASE): return True
    if re.search(r'B[sl][lr]\d', func, re.IGNORECASE): return True
    if re.search(r'^U\d', func, re.IGNORECASE): return True
    if re.search(r'^orf[^_]', func, re.IGNORECASE): return True
    if re.search(r'uncharacterized', func, re.IGNORECASE): return True
    if re.search(r'pseudogene', func, re.IGNORECASE): return True
    if re.search(r'^predicted', func, re.IGNORECASE): return True
    if re.search(r'AGR_', func, re.IGNORECASE): return True
    if re.search(r'similar to', func, re.IGNORECASE): return True
    if re.search(r'similarity', func, re.IGNORECASE): return True
    if re.search(r'glimmer', func, re.IGNORECASE): return True
    if re.search(r'unknown', func, re.IGNORECASE): return True
    if re.search(r'domain', func, re.IGNORECASE): return True
    if re.search(r'^y[a-z]{2,4}\b', func, re.IGNORECASE): return True
    if re.search(r'complete', func, re.IGNORECASE): return True
    if re.search(r'ensang', func, re.IGNORECASE): return True
    if re.search(r'unnamed', func, re.IGNORECASE): return True
    if re.search(r'EG:', func, re.IGNORECASE): return True
    if re.search(r'orf\d+', func, re.IGNORECASE): return True
    if re.search(r'RIKEN', func, re.IGNORECASE): return True
    if re.search(r'Expressed', func, re.IGNORECASE): return True
    if re.search(r'[a-zA-Z]{2,3}\|', func, re.IGNORECASE): return True
    if re.search(r'predicted by Psort', func, re.IGNORECASE): return True
    if re.search(r'^bh\d+', func, re.IGNORECASE): return True
    if re.search(r'cds_', func, re.IGNORECASE): return True
    if re.search(r'^[a-z]{2,3}\d+[^:\+\-0-9]', func, re.IGNORECASE): return True
    if re.search(r'similar to', func, re.IGNORECASE): return True
    if re.search(r' identi', func, re.IGNORECASE): return True
    if re.search(r'ortholog of', func, re.IGNORECASE): return True
    if re.search(r'ortholog of', func, re.IGNORECASE): return True
    if re.search(r'structural feature', func, re.IGNORECASE): return True
    if re.search(r'Phage protein', func, re.IGNORECASE): return True
    if re.search(r'mobile element', func, re.IGNORECASE): return True

    return False

def count_hypothetical_proteins(gbk_file):
    count = 0
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if "product" in feature.qualifiers:
                    # Take the first entry of the 'product' list
                    fn = feature.qualifiers["product"][0]
                    if is_hypothetical(fn):
                        count += 1
    return count

def iterate_genbank_files(genbank_directory, output_csv):
    hypothetical_counts_dict = {}

    for file_name in os.listdir(genbank_directory):
        if file_name.endswith(".gbk") or file_name.endswith(".gb"):
            file_path = os.path.join(genbank_directory, file_name)
            hypothetical_count = count_hypothetical_proteins(file_path)
            hypothetical_counts_dict[file_name] = hypothetical_count

    # Write to CSV
    with open(output_csv, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Filename', 'Hypothetical_Protein_Count'])
        
        for file_name, hypothetical_count in hypothetical_counts_dict.items():
            writer.writerow([file_name, hypothetical_count])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script that takes a directory containing genbank files and writes gene presence absence table")
    parser.add_argument('-d', '--directory', dest='directory', help='Enter the directory containing the genbank files')
    parser.add_argument('-o', dest='output', help='Enter the output tabular format')
    results = parser.parse_args()
    iterate_genbank_files(results.directory, results.output)
