"""
#taken from https://raw.githubusercontent.com/linsalrob/EdwardsLab/
Given a directory of fasta files, do some magic to figure out which ones should be reverse complements
"""

import os
import sys
import argparse

from roblib import bcolors
import re
from roblib import stream_fasta, rc
from itertools import product
from scipy.spatial import distance
import pandas as pd

def count_kmers(dir, kmersize, outdir, verbose=False):
    """
    Count the kmers in all the fasta files
    :param dir: the directory
    :param kmersize: the kmer size
    :param outdir: the output directory
    :param verbose:
    :return:
    """

    bases = {'A', 'C', 'T', 'G'}
    fcount = {}
    rcount = {}

    kmers = ["".join(s) for s in product(bases, repeat=kmersize)]
    fare = re.compile('.fasta$|.fna$|.fa$')
    for f in filter(fare.search, os.listdir(dir)):
        if args.v:
            sys.stderr.write(f"{bcolors.GREEN}Reading {f}{bcolors.ENDC}\n")
        for seqid, seq in stream_fasta(os.path.join(dir, f)):
            fcount[seqid] = []
            rcount[seqid] = []
            for k in kmers:
                fcount[seqid].append(seq.count(k))
                rcount[seqid].append(seq.count(rc(k)))

    # now we have a hash with all the fwd counts and fwd->rev counts. Which is closer
    dists = {}
    for s in fcount:
        dists[s] = {}
        for t in fcount:
            distf = distance.euclidean(fcount[s], fcount[t])
            distr = distance.euclidean(fcount[s], rcount[t])
            deltad = distf/(distr+0.000001) # add episilon incase rd =  0

            """
            # Uncomment this line to print the results as a list
            print("\t".join(map(str, [s, t, fd, rd, dd])))
            """
            """
            deltad is the ratio between the number of fwd kmers and the 
            number of reverse kmers. When deltad == 0, the sequences
            are the same.
            When deltad is 1 the sequences are random
            When deltad is very high (>2 but sometimes much bigger) the 
            sequences are the reverse complement 
            
            We set:
            -1 -- they are the same
             0 -- no relationship
             1 -- need reverse complementing
            """

            if deltad < 0.5:
                dists[s][t] = -1
            elif deltad  > 2:
                dists[s][t] = 1
            else:
                dists[s][t] = 0

    # now we just need to figure out which ones minimize the score by reverse complementing them
    # we are going to do this with pandas. I have a jupyter notebook (calculate_which_to_rc) about
    # this
    df = pd.DataFrame.from_dict(dists)
    seen = {}
    oldsum = 0
    maxiters = 100
    iters = 0

    """
    there are three ways to control this loop:
     - run it for n (100) times
     - run it until we see the same thing again
     - run it while we are always improving.
     
    It probably doesn't matter which one, we just run it 100 times
    We also check whether we are just flipping the same thing
    """
    oldsum = 0
    seen = {}
    maxiters = 100
    iters = 0
    lastflip = None
    while True:
        # this is based on number of runs
        iters += 1
        if iters > maxiters:
            break

        totalsum = sum(df.sum(axis=1))
        """
        # uncomment to run while we are getting better
        if totalsum > oldsum:
            sys.stderr.write(f"Total sum is {totalsum}. Oldsum is {oldsum}. Breaking\n")
            break
        oldsum = totalsum
        """

        ## what is the top one to change
        idx = df.sum(axis=1).sort_values(ascending=False).head(2).index
        nm = idx[0]
        if nm == lastflip and len(idx) > 1:
            nm = idx[1]
        lastflip = nm

        """
        # uncomment to break on if we have run before
        if nm in seen:
            sys.stderr.write(f"Total sum is {totalsum}. Found {nm} ({seen[nm]}) again. Breaking\n")
            break
        """

        seen[nm] = seen.get(nm, 0) + 1
        df.loc[nm] *= -1
        df[nm] *= -1
    # finally, we need to know which sequences to reverse complement. Only those where seen is an odd number!

    torc = set()
    for s in seen:
        if seen[s] % 2:
            torc.add(s)

    if args.v:
        sys.stderr.write(f"{bcolors.BLUE}Will reverse complement: {torc}{bcolors.ENDC}\n")
    # and now we read all the fasta files again and put them in output directory
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for f in filter(fare.search, os.listdir(dir)):
        if args.v:
            sys.stderr.write(f"{bcolors.GREEN}Re-reading to check for rc {f}{bcolors.ENDC}\n")
        with open(os.path.join(outdir, f), 'w') as out:
            for seqid, seq in stream_fasta(os.path.join(dir, f)):
                if seqid in torc:
                    out.write(f">{seqid}_rc\n{rc(seq)}\n")
                else:
                    out.write(f">{seqid}\n{seq}\n")




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot a heatmap")
    parser.add_argument('-d', help='directory of fasta files', required=True)
    parser.add_argument('-k', help='kmer size', required=True, type=int)
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    count_kmers(args.d, args.k, args.o, args.v)