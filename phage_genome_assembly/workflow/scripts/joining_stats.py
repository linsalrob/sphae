#!/usr/bin/env python

import sys 
import os
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from itertools import chain

def consolidate(coverage, viral_check, graph_comp, checkv, out):
    cov=pd.read_csv(coverage, sep="\t")
    vv=pd.read_csv(viral_check)
    gcp=pd.read_csv(graph_comp, sep="\t")
    cv=pd.read_csv(checkv, sep="\t")
    tmp=cov.join(gcp.set_index("ContigID"), on ="Contig")
    tmp2=pd.merge(tmp, vv, how="inner", left_on="Contig", right_on="Contig name")
    o=pd.merge(tmp2, cv, how="inner", left_on="Contig", right_on="contig_id")
    #exporting to file, csv
    o.to_csv(out, encoding='utf-8')


consolidate(snakemake.input.coverm, snakemake.input.viralverify, snakemake.input.comp, snakemake.input.checkv, snakemake.output.csv)





# if __name__=='__main__' :
#     parser=argparse.ArgumentParser(description="Joining the stats tables with converage, graph components and viralverify results")
#     parser.add_argument ('-c', dest='coverm', help='Enter the coverm result filename')
#     parser.add_argument ('-v', dest='viralverify', help= 'Enter the viralverify results file')
#     parser.add_argument ('-g', dest='graph', help= 'Enter the graph components tsv file')
#     parser.add_argument ('-o', dest='output', help= 'Enter the output file name')
#     results=parser.parse_args()
#     consolidate(results.coverm, results.viralverify, results.graph, results.output)