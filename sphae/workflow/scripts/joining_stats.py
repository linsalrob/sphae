#!/usr/bin/env python

import sys 
import os
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from itertools import chain

def consolidate(viral_check, graph_comp, checkv, out):
    # Check file sizes
    file_sizes = {
        "viral_check": os.path.getsize(viral_check),
        "graph_comp": os.path.getsize(graph_comp),
        "checkv": os.path.getsize(checkv)
    }

    # If any input file is empty, touch the output file and return
    if any(size == 0 for size in file_sizes.values()):
        open(out, 'a').close()
        return
    
    vv=pd.read_csv(viral_check)
    gcp=pd.read_csv(graph_comp, sep="\t")
    cv=pd.read_csv(checkv, sep="\t")
    tmp=pd.merge(gcp, vv, left_on="ContigID", right_on="Contig name")
    o=pd.merge(tmp, cv, how="inner", left_on="ContigID", right_on="contig_id")
    #exporting to file, csv
    o.to_csv(out, encoding='utf-8')


consolidate(snakemake.input.viralverify, snakemake.input.comp, snakemake.input.checkv, snakemake.output.csv)

#if __name__=='__main__' :
#    parser=argparse.ArgumentParser(description="Joining the stats tables with converage, graph components and viralverify results")
#    parser.add_argument ('-v', dest='viralverify', help= 'Enter the viralverify results file')
#    parser.add_argument ('-g', dest='graph', help= 'Enter the graph components tsv file')
#    parser.add_argument ('-c', dest='checkv', help='Enter the checkv results')
#    parser.add_argument ('-o', dest='output', help= 'Enter the output file name')
#    results=parser.parse_args()
#    consolidate(results.viralverify, results.graph, results.checkv, results.output)