#!/usr/bin/env python

import sys
import numpy as np
import re
import os
import pandas as pd
from os.path import exists
from collections import defaultdict
import argparse

# Define the picking_contigs function
def picking_contigs(file,out):
    data = pd.DataFrame()
    if (os.path.exists(file) and os.path.getsize(file) > 0):
        data = pd.read_csv(file, header=0)
        datav = data[data["Length_x"] > 1000]
        datav = datav[datav["Prediction"] == "Virus"]
        datac = datav[datav["completeness"]> 70.00]
        #print (len(data))
        #print (data)
    else:
        open(out, 'a').close()
        return None

    if (len(datac))==0:
        print("Genome wasn't assembled well")
        datac.to_csv(out, index=False)
        #return None
            
    elif (len(datac))>1:
        if (datac["Connections"] > 0).any():
            print ("The genome is fragmented")
            datac.to_csv(out, index=False)
        #return None
    
    elif (len(datac))==1:
        #print ("entering this if statement")
        if (datac["Connections"] == 0).any():
            print("The genome is assembled, yay!")
            datac.to_csv(out, index=False)
        else:
            print ("The genome has some regions that are frgamented, but mostly assembled")
            print ("Take a look at the assembly graph file in bandage for more information on where the genome is fragmented")
            datac.to_csv(out, index=False)


picking_contigs(snakemake.input.csv, snakemake.output.out)

#if __name__=='__main__' :
#    parser=argparse.ArgumentParser(description="Picking the contig candidates from the resulting stats file ")
#    parser.add_argument ('-c', dest='file', help='Enter the stats result filename')
#    parser.add_argument ('-o', dest='out', help= 'Enter the output file name')
#    results=parser.parse_args()
#    picking_contigs(results.file, results.out)