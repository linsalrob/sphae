#!/usr/bin/env python

import sys 
import os
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from itertools import chain

def consolidate(coverage, viral_check, graph_comp, out):
    cov=pd.read_csv(coverage, sep="\t")
    vv=pd.read_csv(viral_check)
    gcp=pd.read_csv(graph_comp, sep="\t")
    tmp=cov.join(gcp.set_index("ContigID"), on ="Contig")
    #print (vv)
    o=pd.merge(tmp, vv, how="inner", left_on="Contig", right_on="Contig name")
    
    #exporting to file, csv
    o.to_csv(out, encoding='utf-8')



if __name__=='__main__' :
    parser=argparse.ArgumentParser(description="Joining the stats tables with converage, graph components and viralverify results")
    parser.add_argument ('-c', dest='coverm', help='Enter the coverm result filename')
    parser.add_argument ('-v', dest='viralverify', help= 'Enter the viralverify results file')
    parser.add_argument ('-g', dest='graph', help= 'Enter the graph components tsv file')
    parser.add_argument ('-o', dest='output', help= 'Enter the output file name')
    results=parser.parse_args()
    consolidate(results.coverm, results.viralverify, results.graph, results.output)