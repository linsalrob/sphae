"""
Parse a blast file and create a blast result object
"""

import os
import sys
import argparse
import gzip

class BlastResult():
    def __init__(self, query, db, percent_id, alignment_length, gaps, mismatches, query_start, query_end, db_start, db_end,
                 evalue, bitscore, query_length=None, subject_length=None):
        self.query = query
        self.db = db
        self.alignment_length = int(alignment_length)
        self.percent_id = float(percent_id)
        self.gaps = int(gaps)
        self.mismatches = int(mismatches)
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.db_start = int(db_start)
        self.db_end = int(db_end)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)
        if query_length:
            self.query_length = int(query_length)
        if subject_length:
            self.subject_length = int(subject_length)

    def is_significant(self):
        """
        Is this hit significant
        :return: boolean
        """

        return self.evalue < 1e-5


def stream_blast_results(blastf, verbose=False):
    """
    Parse a tab-separated blast file and stream the results
    :param blastf: the file to stream
    :return: a stream of BlastResults
    """

    if blastf.endswith('.gz'):
        qin = gzip.open(blastf, 'rt')
    else:
        qin = open(blastf, 'r')

    while True:
        l = qin.readline()
        if not l:
            break
        p = l.strip().split("\t")
        #if verbose:
        #    sys.stderr.write("LINE: {}\n".format(p))
        br = BlastResult(*p)
        yield br

    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', help='')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()
