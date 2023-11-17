
import sys

"""
Translate and back translate DNA to protein
"""


aa_1_to_3_letter ={"A" : "Ala", "C" : "Cys", "D" : "Asp", "E" : "Glu", "F" : "Phe", "G" : "Gly", "H" : "His",
                   "I" : "Ile", "K" : "Lys", "L" : "Leu", "M" : "Met", "N" : "Asn", "P" : "Pro", "Q" : "Gln",
                   "R" : "Arg", "S" : "Ser", "T" : "Thr", "V" : "Val", "W" : "Trp", "Y" : "Tyr", "*" : "Stop"}

aa_1_letter_order = "A C D E F G H I K L M N P Q R S T V W Y".split()  # Alpha by 1 letter
aa_3_letter_order = "A R N D C Q E G H I L K M F P S T W Y V".split()  # PAM matrix order
aa_n_codon_order  = "L R S A G P T V I C D E F H K N Q Y M W".split()

genetic_code = {
    # DNA version
    "TTT": 'F',  "TCT": 'S',  "TAT": 'Y',  "TGT": 'C',
    "TTC": 'F',  "TCC": 'S',  "TAC": 'Y',  "TGC": 'C',
    "TTA": 'L',  "TCA": 'S',  "TAA": '*',  "TGA": '*',
    "TTG": 'L',  "TCG": 'S',  "TAG": '*',  "TGG": 'W',
    "CTT": 'L',  "CCT": 'P',  "CAT": 'H',  "CGT": 'R',
    "CTC": 'L',  "CCC": 'P',  "CAC": 'H',  "CGC": 'R',
    "CTA": 'L',  "CCA": 'P',  "CAA": 'Q',  "CGA": 'R',
    "CTG": 'L',  "CCG": 'P',  "CAG": 'Q',  "CGG": 'R',
    "ATT": 'I',  "ACT": 'T',  "AAT": 'N',  "AGT": 'S',
    "ATC": 'I',  "ACC": 'T',  "AAC": 'N',  "AGC": 'S',
    "ATA": 'I',  "ACA": 'T',  "AAA": 'K',  "AGA": 'R',
    "ATG": 'M',  "ACG": 'T',  "AAG": 'K',  "AGG": 'R',
    "GTT": 'V',  "GCT": 'A',  "GAT": 'D',  "GGT": 'G',
    "GTC": 'V',  "GCC": 'A',  "GAC": 'D',  "GGC": 'G',
    "GTA": 'V',  "GCA": 'A',  "GAA": 'E',  "GGA": 'G',
    "GTG": 'V',  "GCG": 'A',  "GAG": 'E',  "GGG": 'G',

    #  The following ambiguous encodings are not necessary,  but
    #  speed up the processing of some ambiguous triplets: 

    "TTY": 'F',  "TCY": 'S',  "TAY": 'Y',  "TGY": 'C',
    "TTR": 'L',  "TCR": 'S',  "TAR": '*',
    "TCN": 'S',
    "CTY": 'L',  "CCY": 'P',  "CAY": 'H',  "CGY": 'R',
    "CTR": 'L',  "CCR": 'P',  "CAR": 'Q',  "CGR": 'R',
    "CTN": 'L',  "CCN": 'P',  "CGN": 'R',
    "ATY": 'I',  "ACY": 'T',  "AAY": 'N',  "AGY": 'S',
    "ACR": 'T',  "AAR": 'K',  "AGR": 'R',
    "ACN": 'T',
    "GTY": 'V',  "GCY": 'A',  "GAY": 'D',  "GGY": 'G',
    "GTR": 'V',  "GCR": 'A',  "GAR": 'E',  "GGR": 'G',
    "GTN": 'V',  "GCN": 'A',  "GGN": 'G'
}

amino_acid_codons_DNA = {
        "L": "TTA TTG CTA CTG CTT CTC".split(),
        "R": "AGA AGG CGA CGG CGT CGC".split(),
        "S": "AGT AGC TCA TCG TCT TCC".split(),
        "A": "GCA GCG GCT GCC".split(),
        "G": "GGA GGG GGT GGC".split(),
        "P": "CCA CCG CCT CCC".split(),
        "T": "ACA ACG ACT ACC".split(),
        "V": "GTA GTG GTT GTC".split(),
        "I": "ATA ATT ATC".split(),
        "C": "TGT TGC".split(),
        "D": "GAT GAC".split(),
        "E": "GAA GAG".split(),
        "F": "TTT TTC".split(),
        "H": "CAT CAC".split(),
        "K": "AAA AAG".split(),
        "N": "AAT AAC".split(),
        "Q": "CAA CAG".split(),
        "Y": "TAT TAC".split(),
        "M": "ATG".split(),
        "U": "TGA".split(),
        "W": "TGG".split(),
        "l": "tta ttg cta ctg ctt ctc".split(),
        "r": "aga agg cga cgg cgt cgc".split(),
        "s": "agt agc tca tcg tct tcc".split(),
        "a": "gca gcg gct gcc".split(),
        "g": "gga ggg ggt ggc".split(),
        "p": "cca ccg cct ccc".split(),
        "t": "aca acg act acc".split(),
        "v": "gta gtg gtt gtc".split(),
        "i": "ata att atc".split(),
        "c": "tgt tgc".split(),
        "d": "gat gac".split(),
        "e": "gaa gag".split(),
        "f": "ttt ttc".split(),
        "h": "cat cac".split(),
        "k": "aaa aag".split(),
        "n": "aat aac".split(),
        "q": "caa cag".split(),
        "y": "tat tac".split(),
        "m": "atg".split(),
        "u": "tga".split(),
        "w": "tgg".split(),
        '*': "TAA TAG TGA".split()
}


def translate_dna(sequence, verbose=False):
    """
    Translate a DNA sequence and return a protein string
    :param sequence: The DNA sequence to translate
    :param verbose: More output
    :return: a protein string
    """

    posn=0
    trans=""
    while posn < len(sequence)-3:
        codon = sequence[posn:posn+3]
        if codon not in genetic_code:
            sys.stderr.write("Uknown codon: {}\n".format(codon))
            trans += "X"
            continue
        trans += genetic_code[codon]
        posn += 3
    return trans

def print_codon_table(ambiguous=False):
    """
    Print the codon usage table
    :param ambiguous: Print codons with ambiguous bases
    :return:
    """

    c = sorted(genetic_code.keys())
    for codon in c:
        test = codon.replace('A', '')
        test = test.replace('G', '')
        test = test.replace('C', '')
        test = test.replace('T', '')
        if test and not ambiguous:
            continue

        print("{}\t{}".format(codon, aa_1_to_3_letter[genetic_code[codon]]))
