#!/usr/bin/env python

"""
code to rename the genomes with the sample name
genome name - sampleName_contig001... 

Also add keep a tsv file keep the old contig name and the new name

Code to run after genome selection prior running annotation
"""

import sys
import re
import os
from Bio import SeqIO
import argparse
import csv

#file - input fasta file; sample - sample name; out - output file name 
def rename_contigs(filein,sample,out,dicts):
    fin = open(filein, "r")
    fout=open(out,"w")
    Records=[]
    count=0
    keep_track={} #keys are the new names and the values are the old names 
    for record in SeqIO.parse(fin,'fasta'):
        oldID=record.id
        record.id = sample + "_contig_" + str(count)
        keep_track[record.id]=oldID
        Records.append(record)
        count=count+1

    for sequence in Records:
        SeqIO.write(sequence, fout,"fasta")

    with open(dicts, 'w') as csvfile:
        for key in keep_track.keys():
            csvfile.write("%s, %s\n" % (key, keep_track[key]))

rename_contigs(snakemake.input.fin, snakemake.params.s, snakemake.output.out, snakemake.output.csv)
#parser = argparse.ArgumentParser()

#if __name__=='__main__' :
#    parser=argparse.ArgumentParser(description="Rename the genomes with the sample name and the contig ID")
#    parser.add_argument ('-i', dest='file', help='Enter the genome fasta file')
#    parser.add_argument ('-s', dest='sample', help='Enter the sample name')
#    parser.add_argument ('-o', dest='out', help='Enter the output fasta file name')
#    parser.add_argument ('-t', dest='csv', help= 'Enter the csv filename to keep track of the name change')
#    results=parser.parse_args()
#    rename_contigs(results.file, results.sample, results.out, results.csv)