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
    if (os.path.exists(file) == True):
        data = pd.read_csv(file, header=0)
        data = data[data["Length_x"] > 1000]
        data = data[data["Prediction"] == "Virus"]
        data = data[data["Mean"] > 1]
        data = data[data["completeness"]> 90.00]
        #print (len(data))
        print (data)

    if (len(data))==0:
        print("Genome wasn't assembled well")
        return None
            
    elif (len(data))>1:
        if (data["Connections"] > 0).any():
            print ("The genome is fragmented")
            print ("Subsample the reads and run again, to do this add the parameter in the config file")
        return None
    
    elif (len(data))==1:
        #print ("entering this if statement")
        if (data["Connections"] == 0).any():
            print("The genome is assembled, yay!")
            data.to_csv(out, index=False)
        else:
            print ("The genome has some regions that are frgamented, but mostly assembled")
            print ("Take a look at the assembly graph file in bandage for more information on where the genome is fragmented")
            data.to_csv(out, index=False)


picking_contigs(snakemake.input.csv, snakemake.output.out)

#if __name__=='__main__' :
#    parser=argparse.ArgumentParser(description="Picking the contig candidates from the resulting stats file ")
#    parser.add_argument ('-c', dest='file', help='Enter the stats result filename')
#    parser.add_argument ('-o', dest='out', help= 'Enter the output file name')
#    results=parser.parse_args()
#    picking_contigs(results.file, results.out)