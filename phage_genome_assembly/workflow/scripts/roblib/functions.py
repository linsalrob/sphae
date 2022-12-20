import re


"""

Determine whether a protein function is hypothetical.

"""

def is_hypothetical(func):
    """
    Returns True if the function is hypothetical. Otherwise returns false
    :param func: string
    :return: boolean
    """

    if not func: return True
    if func.lower() == 'hypothetical protein': return True
    if re.search('lmo\d+ protein', func, re.IGNORECASE): return True
    if re.search('hypoth', func, re.IGNORECASE): return True
    if re.search('conserved protein', func, re.IGNORECASE): return True
    if re.search('gene product', func, re.IGNORECASE): return True
    if re.search('interpro', func, re.IGNORECASE): return True
    if re.search('B[sl][lr]\d', func, re.IGNORECASE): return True
    if re.search('^U\d', func, re.IGNORECASE): return True
    if re.search('^orf[^_]', func, re.IGNORECASE): return True
    if re.search('uncharacterized', func, re.IGNORECASE): return True
    if re.search('pseudogene', func, re.IGNORECASE): return True
    if re.search('^predicted', func, re.IGNORECASE): return True
    if re.search('AGR_', func, re.IGNORECASE): return True
    if re.search('similar to', func, re.IGNORECASE): return True
    if re.search('similarity', func, re.IGNORECASE): return True
    if re.search('glimmer', func, re.IGNORECASE): return True
    if re.search('unknown', func, re.IGNORECASE): return True
    if re.search('domain', func, re.IGNORECASE): return True
    if re.search('^y[a-z]{2,4}\b', func, re.IGNORECASE): return True
    if re.search('complete', func, re.IGNORECASE): return True
    if re.search('ensang', func, re.IGNORECASE): return True
    if re.search('unnamed', func, re.IGNORECASE): return True
    if re.search('EG:', func, re.IGNORECASE): return True
    if re.search('orf\d+', func, re.IGNORECASE): return True
    if re.search('RIKEN', func, re.IGNORECASE): return True
    if re.search('Expressed', func, re.IGNORECASE): return True
    if re.search('[a-zA-Z]{2,3}\|', func, re.IGNORECASE): return True
    if re.search('predicted by Psort', func, re.IGNORECASE): return True
    if re.search('^bh\d+', func, re.IGNORECASE): return True
    if re.search('cds_', func, re.IGNORECASE): return True
    if re.search('^[a-z]{2,3}\d+[^:\+\-0-9]', func, re.IGNORECASE): return True
    if re.search('similar to', func, re.IGNORECASE): return True
    if re.search(' identi', func, re.IGNORECASE): return True
    if re.search('ortholog of', func, re.IGNORECASE): return True
    if re.search('ortholog of', func, re.IGNORECASE): return True
    if re.search('structural feature', func, re.IGNORECASE): return True
    if re.search('Phage protein', func, re.IGNORECASE): return True
    if re.search('mobile element', func, re.IGNORECASE): return True



    return False





if __name__ == "__main__":

    for fn in ['Real Function', 'lmo24 protein', 'Hypothetical Protein']:
        print("{}\t{}".format(fn, is_hypothetical(fn)))