import os
import sys
import string
__author__ = 'Rob Edwards'




def ascii_clean(s):
    """Remove non-ascii characters from a string"""
    return filter(lambda x: x in string.printable, s)

