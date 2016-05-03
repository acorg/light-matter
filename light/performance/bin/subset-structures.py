#!/usr/bin/env python

"""
Read a fasta file from stdin, and write all possible substrings to stdout.
Note that the range of lengths of the new substrings is defined between the
shortest sequence and the longest sequence in the file.
"""

import sys

from dark.fasta import FastaReads
from dark.reads import AAReadWithX


def makeSubstring(length, string):
    """
    Take a string and make substrings of a specified length.

    @param length: An C{int} of length of the substrings.
    @param string: A {dark.AAReadsWithX} instance of a sequence to make
        substrings from.
    """
    assert(length < len(string))
    numberOfStrings = len(string) - length + 1
    for i in range(numberOfStrings):
        newStringId = '%s[%d:%d]' % (string.id, i, length + i)
        yield AAReadWithX(newStringId, string.sequence[i:length + i])


allStrings = list(FastaReads(sys.stdin, readClass=AAReadWithX,
                             checkAlphabet=0))

allLengths = [len(string) for string in allStrings]
minLength = min(allLengths)
maxLength = max(allLengths)

for length in range(minLength, maxLength + 1):
    for string in allStrings:
        if length == len(string):
            print(string.toString(format_='fasta'), end='')
        elif length < len(string):
            for newString in makeSubstring(length, string):
                print(newString.toString(format_='fasta'), end='')
