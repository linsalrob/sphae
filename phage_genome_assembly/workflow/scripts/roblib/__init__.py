import os
import sys

__author__ = 'Rob Edwards'
from .stats import mean, median, stdev
from .sequences import read_fasta, readFasta, stream_fastq, stream_fasta, stream_paired_fastq, stream_gfa_sequences
from .sequences import write_fastq, qual_to_numbers
from .dna import rc, shannon
from .geography import latlon2distance
from .strings import ascii_clean
from .functions import is_hypothetical
from .newick import Newick_Tree
from .dnadist import parse_dnadist
from .blast import stream_blast_results
from .translate import translate_dna
from .bcolors import bcolors
from .rob_error import SequencePairError, FastqFormatError
from .colours import colours, colors, message
from .genbank import genbank_to_faa, genbank_to_fna, genbank_to_orfs, genbank_seqio
from .genbank import genbank_to_ptt, genbank_to_functions, feature_id, genbank_to_pandas
from .file_chooser import choose_a_file, write_a_file

__all__ = [
    'mean', 'median', 'stdev',
    'read_fasta', 'readFasta', 'stream_fastq', 'stream_fasta', 'stream_paired_fastq', 'stream_gfa_sequences',
    'write_fastq', 'qual_to_numbers',
    'rc', 'shannon',
    'latlon2distance',
    'ascii_clean', 'is_hypothetical', 'Newick_Tree', 'parse_dnadist',
    'stream_blast_results',
    'translate_dna',
    'bcolors', 'colours', 'colors', 'message',
    'SequencePairError', 'FastqFormatError',
    'genbank_to_faa', 'genbank_to_fna', 'genbank_to_orfs', 'genbank_to_ptt', 'genbank_seqio', 'genbank_to_functions',
    'feature_id', 'genbank_to_pandas',
    'choose_a_file', 'write_a_file'
    ]
