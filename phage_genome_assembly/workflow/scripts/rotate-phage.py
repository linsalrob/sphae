"""
Taken from https://raw.githubusercontent.com
Rotate a fasta sequence so that the start is just before the terminase.

We need a blast hit to the terminase large subunit, but we don't necessarily need
a hit to the terminase small subunit. These two are usually adjacent and on the same
strand (in 2,127 out of 2,165 times!), and so we look upstream of the
large subunit. If we have the small subunit and it is upstream of the
large subunit, we will put the break before the small subunit. Otherwise
we put the break before the large subunit.
"""

import os
import sys
import argparse
from roblib import stream_blast_results, rc, read_fasta
from roblib import stream_fasta, rc

def read_blast_file(bf, evalue, verbose=False):
    """
    Read the blast file and store tuples of contig, start, stop for each hit
    :param bf: blast file to parse
    :param evalue: E value to keep matches
    :param verbose: more output
    :return:
    """

    results = {}

    for br in stream_blast_results(bf, verbose=verbose):
        if br.evalue > evalue:
            continue
        if br.query not in results:
            results[br.query] = { '-ve' : set(), '+ve' : set() }
        strand = '+ve'
        (s, e) = (br.query_start, br.query_end)
        if br.query_end < br.query_start:
            strand = '-ve'
            (s, e) = (br.query_end, br.query_start)
            
        results[br.query][strand] = results[br.query][strand] | set(range(s, e))
    return results


def is_consecutive(l):
    """
    Determines whether a list contains consecutive numbers. If it does, the sum of the
    numbers should be the same as n * (min + max) / 2 (where n is num elements)
    :param l: A set of numbers
    :return: True or False
    """

    total = sum(l)
    mathsum = len(l) * (min(l) + max(l)) / 2
    return total == mathsum

def find_gene(brs, force, verbose=False):
    """
    Test for where the genes are
    :param brs: The blast results
    :param force: force a gene, even if there are multiple discontinuous genes
    :param verbose: more output
    :return:
    """

    genes = {}
    contigs = brs.keys()
    if len(contigs) == 1:
        sys.stderr.write("There was only one contig with a match ({})\n".format(list(contigs)[0]))
    else:
        sys.stderr.write("Multiple contigs had matches: {}\n".format(contigs))

    for c in contigs:
        sys.stderr.write(f"{c}\n")
        for strand in ['-ve', '+ve']:
            if len(brs[c][strand]) == 0:
                continue
            if is_consecutive(brs[c][strand]):
                sys.stderr.write("There is a consecutive match to {} {} {}\n".format(c, max(brs[c][strand]), min(brs[c][strand])))
                if c not in genes:
                    genes[c] = {strand : [min(brs[c][strand]), max(brs[c][strand])] }
                else:
                    sys.stderr.write("Warning: two hits on different strands for {}\n".format(c))
                    genes[c][strand]  = [min(brs[c][strand]), max(brs[c][strand])]
            else:
                if verbose:
                    startstops = []
                    start = min(brs[c][strand])
                    stop = start
                    curr = start
                    while curr <= max(brs[c][strand]):
                        #sys.stderr.write("{}\t{}\t{}\n".format(curr, start, stop, curr in brs[c][strand]))
                        if curr in brs[c][strand]:
                            stop = curr
                            if not start:
                                start = curr
                        else:
                            if start:
                                startstops.append([start, stop])
                                start = None
                                stop = None
                        curr += 1
                    if start:
                        startstops.append([start, stop])

                    sys.stderr.write("Min: {} Max: {}\n".format(min(brs[c][strand]), max(brs[c][strand])))
                    sys.stderr.write("List of starts and stops: {}\n".format(startstops))
                sys.stderr.write("There are multiple discontinuous matches to {}. Try adjusting the evalue\n".format(c))
                if force:
                    sys.stderr.write("***** WARNING: We are forcing a gene from {} to {}\n".format(min(brs[c][strand]), max(brs[c][strand])))
                    if c not in genes:
                        genes[c] = {strand: [min(brs[c][strand]), max(brs[c][strand])]}
                    else:
                        sys.stderr.write("Warning: two hits on different strands for {}\n".format(c))
                        genes[c][strand] = [min(brs[c][strand]), max(brs[c][strand])]

    return genes


def introduce_break(fastaf, lsgenes, ssgenes, numbp, verbose=False):
    """
    Introduce a break into the contigs upstream of the large subunit, but also ignore the small subunit
    :param fastaf: fasta file to break
    :param lsgenes: dict of ls genes
    :param ssgenes: dict of ss genes
    :param numbp: how far upstream to break the DNA
    :param verbose: more output
    :return:
    """

    for seqdef, seq in stream_fasta(fastaf):
        seqid = seqdef.split(" ")[0]
        if verbose:
            sys.stderr.write(f"Looking for {seqid}\n")
        if seqid in lsgenes:
            breakat = 0
            dorc = False
            if '+ve' in lsgenes[seqid]:
                # the terminase is on the + strand so we go n bp less than the mininmum
                breakat = lsgenes[seqid]['+ve'][0] - numbp
                # check to see if the small subunit is here too
                if seqid in ssgenes and '+ve' in ssgenes[seqid]:
                    if ssgenes[seqid]['+ve'][0] - numbp < breakat:
                        breakat = ssgenes[seqid]['+ve'][0] - numbp
                        if verbose:
                            sys.stderr.write("Because we have a small subunit at {} we moved the start\n".format(ssgenes[seqid]['+ve'][0]))
            elif '-ve' in lsgenes[seqid]:
                dorc = True
                # the terminase is on the - strand so we go n bp more than the maximum
                breakat = lsgenes[seqid]['-ve'][1] + numbp
                # check to see if the small subunit is here too
                if seqid in ssgenes and '-ve' in ssgenes[seqid]:
                    if ssgenes[seqid]['-ve'][1] + numbp > breakat:
                        breakat = ssgenes[seqid]['-ve'][1] + numbp
                        if verbose:
                            sys.stderr.write("Because we have a small subunit at {} we moved the start\n".format(ssgenes[seqid]['-ve'][0]))
            else:
                sys.stderr.write("No terminase found in {}\n".format(seqid))

            newseq = seq
            if breakat > 0:
                sys.stderr.write("Reformating the sequence {}\n".format(seqid))
                newseq = seq[breakat:] + seq[0:breakat]
                if dorc:
                    newseq = rc(newseq)
            print(">{}\n{}".format(seqdef, newseq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Recircularize the genome at the terminase")
    parser.add_argument('-f', help='fasta file of the genome', required=True)
    parser.add_argument('-l', help='blast output compared to terminase LARGE gene', required=True)
    parser.add_argument('-s', help='blast output compared to terminase SMALL gene')
    parser.add_argument('-e', help='evalue cutoff (default=all)', type=int, default=100)
    parser.add_argument('-n', help='number of bp from the end of the terminase to introduce the break (default=100)', type=int, default=100)
    parser.add_argument('--force', help='Force the rotation, even if there is more than one gene. In this case, we use the longest possible gene!', action='store_true')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    ssgenes = {}
    if args.s:
        sys.stderr.write("Parsing small subunit\n")
        ss = read_blast_file(args.s, args.e, args.v)
        ssgenes = find_gene(ss, args.force, args.v)

    sys.stderr.write("Parsing large subunit\n")
    ls = read_blast_file(args.l, args.e, args.v)
    lsgenes = find_gene(ls, args.force,  args.v)

    introduce_break(args.f, lsgenes, ssgenes, args.n, args.v)