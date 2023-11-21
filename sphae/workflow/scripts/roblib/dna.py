import os
import sys

__author__ = 'Rob Edwards'


def rc(dna):
    """
    Reverse complement a DNA sequence

    :param dna: The DNA sequence
    :type dna: str
    :return: The reverse complement of the DNA sequence
    :rtype: str
    """
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq


def shannon(dna, word):
    """
    Count the Shannon entropy for a DNA sequence

    :param dna: The DNA sequence
    :type dna: str
    :param word: Word length
    :type word: int
    :return: the Shannon entropy
    :rtype: float
    """
    count = {}
    for i in range(len(dna) - word):
        substr = dna[i:i + word]
        count[substr] = count.get(substr, 0) + 1

    hp = 0
    p = len(dna) - word
    for s in count:
        n = 1.0 * count[s]
        hp += float((n / p) * log(n / p))
    return 0 - hp

def kmers(dna, k):
    """
    Given a dna sequence return a hash with all kmers of length k and their frequency.

    This method does NOT use the reverse complement, it only checks the strand you supply.

    :param dna: the dna sequence
    :param k: the length of the kmer
    :return: a hash of kmers and abundance
    """

    counts = {}
    for i in range(len(dna)-k):
        subseq=dna[i:i+k]
        counts[subseq] = counts.get(subseq, 0) + 1

    return counts


