import os
import sys
import gzip

import subprocess
from .rob_error import SequencePairError, FastqFormatError
from .colours import colours, message

__author__ = 'Rob Edwards'


def read_fasta(fname: str, whole_id: bool = True, qual: bool = False) -> dict:
    """
    Read a fasta file and return a hash.

    If wholeId is set to false only the first part of the ID
    (upto the first white space) is returned

    :param fname: The file name to read
    :param whole_id: Whether to keep the whole id, or trim to first whitespace (default = all)
    :param qual: these are quality scores (so add a space between lines!)
    :return: dict
    """

    try:
        if fname.endswith('.gz'):
            f = gzip.open(fname, 'rt')
        elif fname.endswith('.lrz'):
            f = subprocess.Popen(['/usr/bin/lrunzip', '-q', '-d', '-f', '-o-', fname], stdout=subprocess.PIPE).stdout
        else:
            f = open(fname, 'r')
    except IOError as e:
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write("Message: \n" + str(e.message) + "\n")
        sys.exit("Unable to open file " + fname)

    seqs = {}
    seq = ""
    seqid = ""
    for line in f:
        line = line.rstrip('\r\n')
        if line.startswith(">"):
            if seqid != "":
                seqs[seqid] = seq
                seq = ""
            seqid = line.replace(">", "", 1)
            if not whole_id and seqid.count(" ") > 0:
                seqids = seqid.split(" ")
                seqid = seqids[0]
        else:
            if qual:
                seq += " " + line
            else:
                seq += line

    seqs[seqid] = seq.strip()
    return seqs


def readFasta(file, whole_id=True):
    """
    Read a fasta file and return a hash.

    If wholeId is set to false only the first part of the ID (upto the first white space) is returned
    """
    return read_fasta(file, whole_id)


def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores.

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """

    if fqfile.endswith('.gz'):
        qin = gzip.open(fqfile, 'rt')
    else:
        qin = open(fqfile, 'r')

    linecounter = 0
    while True:
        header = qin.readline()
        linecounter += 1
        if not header:
            break
        if not header.startswith("@"):
            raise FastqFormatError(f"The file does not appear to be a four-line fastq file at line {linecounter}")
        header = header.strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seqid = seqid.replace('@', '')
        seq = qin.readline().strip()
        linecounter += 1
        qualheader = qin.readline()
        if not qualheader.startswith("+"):
            raise FastqFormatError(f"The file does not appear to be a four-line fastq file at line {linecounter}")
        linecounter += 1
        qualscores = qin.readline().strip()
        linecounter += 1
        header = header.replace('@', '', 1)
        if len(qualscores) != len(seq):
            raise FastqFormatError(f"The sequence and qual scores are not the same length at line {linecounter}")
        yield seqid, header, seq, qualscores



def stream_paired_fastq(fqfile1, fqfile2):
    """Read paired fastq files and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores for the left and right pairs

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.

    Should accomodate both /1 /2 // _1 _2 and [space]1 [space]2
    """

    if fqfile1.endswith('.gz'):
        qin1 = gzip.open(fqfile1, 'rt')
    else:
        qin1 = open(fqfile1, 'r')

    if fqfile2.endswith('.gz'):
        qin2 = gzip.open(fqfile2, 'rt')
    else:
        qin2 = open(fqfile2, 'r')

    linecounter = 0
    while True:
        linecounter += 1
        header1 = qin1.readline().strip()
        header2 = qin2.readline().strip()
        if not header1 and not header2:
            break

        if not header1.startswith("@"):
            raise FastqFormatError(f"The file {fqfile1} does not appear to be a four-line fastq file at line {linecounter}")
        if not header2.startswith("@"):
            raise FastqFormatError(f"The file {fqfile2} does not appear to be a four-line fastq file at line {linecounter}")

        header1 = header1.replace('@', '', 1)
        header2 = header2.replace('@', '', 1)

        seqidparts1 = header1.split(' ')
        seqidparts2 = header2.split(' ')

        seqid1 = seqidparts1[0]
        seqid2 = seqidparts2[0]

        # test for appropriate matching names
        if seqid1.endswith('/1') and seqid2.endswith('/2'):
            seqid1 = seqid1.replace('/1', '')
            seqid2 = seqid2.replace('/2', '')
        elif seqid1.endswith('_1') and seqid2.endswith('_2'):
            seqid1 = seqid1.replace('_1', '')
            seqid2 = seqid2.replace('_2', '')
        elif seqidparts1[1].startswith('1') and seqidparts2[1].startswith('2'):
            # the sequence names contain [space]1 and [space]2
            True
        else:
            raise SequencePairError(f"{colours.RED}We can not match the forward/reverse reads in\n{header1} (seqid: {seqid1})\n{header2} (seqid: {seqid2}){colours.ENDC}\n")

        if seqid1 != seqid2:
            raise SequencePairError(f"{colours.RED}The sequence IDs {seqid1} and {seqid2} do not match and are not paired{colours.ENDC}\n")

        seqid = seqid1.replace('@', '')
        seq1 = qin1.readline().strip()
        linecounter += 1
        qualheader1 = qin1.readline()
        linecounter += 1
        qualscores1 = qin1.readline().strip()
        linecounter += 1


        seq2 = qin2.readline().strip()
        qualheader2 = qin2.readline()
        qualscores2 = qin2.readline().strip()

        if not qualheader1.startswith("+"):
            raise FastqFormatError(f"The file {fqfile1} does not appear to be a four-line fastq file at line {linecounter}")
        if not qualheader2.startswith("+"):
            raise FastqFormatError(f"The file {fqfile2} does not appear to be a four-line fastq file at line {linecounter}")

        if len(qualscores1) != len(seq1):
            raise FastqFormatError(f"The sequence and qual scores in {fqfile1} are not the same length at line {linecounter}")
        if len(qualscores2) != len(seq2):
            raise FastqFormatError(f"The sequence and qual scores in {fqfile2} are not the same length at line {linecounter}")

        yield seqid, header1, seq1, qualscores1, header2, seq2, qualscores2





def stream_fasta(fastafile, whole_id=True):
    """
    Stream a fasta file, one read at a time. Saves memory!

    :param fastafile: The fasta file to stream
    :type fastafile: str
    :param whole_id: Whether to return the whole id (default) or just up to the first white space
    :type whole_id:bool
    :return:A single read
    :rtype:str, str
    """

    try:
        if fastafile.endswith('.gz'):
            f = gzip.open(fastafile, 'rt')
        elif fastafile.endswith('.lrz'):
            f = subprocess.Popen(['/usr/bin/lrunzip', '-q', '-d', '-f', '-o-', fastafile], stdout=subprocess.PIPE).stdout
        else:
            f = open(fastafile, 'r')
    except IOError as e:
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write("Message: \n" + str(e.message) + "\n")
        sys.exit("Unable to open file " + fastafile)

    posn = 0
    while f:
        # first line should start with >
        idline = f.readline()
        if not idline:
            break
        if not idline.startswith('>'):
            sys.exit("Do not have a fasta file at: {}".format(idline))
        if not whole_id:
            idline = idline.split(" ")[0]
        idline = idline.strip().replace('>', '', 1)
        posn = f.tell()
        line = f.readline()
        seq = ""
        while not line.startswith('>'):
            seq += line.strip()
            posn = f.tell()
            line = f.readline()
            if not line:
                break
        f.seek(posn)
        yield idline, seq


def stream_gfa_sequences(gfafile):
    """
    Stream the sequences from a GFA file. At the moment we ignore the 
    rest of the information. This is not supposed to be a parser
    :param gfafile: the gfa file to read
    """

    try:
        if gfafile.endswith('.gz'):
            f = gzip.open(gfafile, 'rt')
        elif gfafile.endswith('.lrz'):
            f = subprocess.Popen(['/usr/bin/lrunzip', '-q', '-d', '-f', '-o-', gfafile], stdout=subprocess.PIPE).stdout
        else:
            f = open(gfafile, 'r')
    except IOError as e:
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write("Message: \n" + str(e.message) + "\n")
        sys.exit("Unable to open file " + gfafile)

    while f:
        l = f.readline()
        if not l:
            break
        if not l.startswith("S"):
            continue
        p = l.strip().split("\t")
        yield p[1], p[2]


def write_fastq(fna, qual, outf, verbose=False):
    """
    Write DNA sequences + quality scores to a fastq file
    :param fna: a dict of the DNA sequences
    :param qual: a dict of the the quality scores
    :param outf: the output file to write
    :param verbose: more output
    :return:
    """

    with open(outf, 'w') as out:
        for k in fna:
            if k not in qual:
                raise FastqFormatError(f"{colours.RED}No quality scores were found for {k}{colours.ENDC}")

            qualstring = None
            if len(fna[k]) == len(qual[k]):
                # the quality is probably already converted
                if verbose:
                    message("Assuming the qaulity scores have already been converted to a string as the seqs and ", "GREEN")
                    message("qual scores are the same length", "GREEN")
                qualstring = qual[k]
            else:
                p = qual[k].strip().split(" ")
                if len(p) != len(fna[k]):
                    msg = f"""
                        {colours.RED}FATAL: For {k} we have sequence:
                        |{fna[k]}|
                        and
                        |{qual[k]}|
                        that became
                        {p}
                        lengths {len(fna[k])} and {len(p)} that are different
                        """
                    raise FastqFormatError(msg)
                qualstring = "".join(map(lambda x: chr(int(x)+33), p))

            out.write(f"@{k}\n{fna[k]}\n+\n{qualstring}\n")

def qual_to_numbers(qualstring):
    """
    Convert the quality string to an array of numbers
    :param qualstring: the quality string
    :return: a list of numbers
    """
    return map(lambda x: ord(x)-33, qualstring)