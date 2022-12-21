"""
Parse the output from Phylip's DNA dist and return the data as a two dimensional array
"""

import os
import sys
import argparse
import re

def parse_dnadist(dnadistancefile):
    """
    Parse a dnadistance file and generate a list of IDs and a two dimensional
    array of the data
    :param dnadistancefile: The dna distance file to parse
    :return: an array of the IDs and a two-dimensional array of the distances
    """

    ids = []
    distarray = []
    distance = {}
    numelements = 0

    data = []

    with open(dnadistancefile, 'r') as f:
        numelements = int(f.readline().strip())
        thisline = []
        for l in f:
            p = re.split("\s+", l.strip())
            if l.startswith(' '):
                thisline.extend(p)
            else:
                if thisline:
                    data.append(thisline)
                thisline = []
                thisline.extend(p)
        if thisline:
            data.append(thisline)

    for p in data:
        sampleid = p.pop(0)
        ids.append(sampleid)
        p = list(map(float, p))
        distarray.append(p)

    if numelements != len(ids):
        sys.stderr.write("ERROR: The dna distance file is supposed to have {} elements but we found {} ids\n".format(numelements, len(ids)))
        sys.exit(-1)

    return ids, distarray


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse the data from a dnadist file')
    parser.add_argument('-d', help='dna distance file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    ids, dists = parse_dnadist(args.d)
    sys.stdout.write("Found {} ids\n".format(len(ids)))
    sys.stdout.write("The distance from the first to second elements is {}\n".format(dists[0][1]))