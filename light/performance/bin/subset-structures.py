#!/usr/bin/env python

"""
Read a fasta file of helices from stdin, and write all possible helix
substrings to stdout.
"""

import sys

from dark.fasta import FastaReads
from dark.reads import AAReadWithX


def makeSubstring(length, helix):
    """
    Take a string and make substrings of a specified length.

    @param length: An C{int} of length of the substrings.
    @param helix: A {dark.AAReadsWithX} instance of sequence to make substrings
        from.
    """
    assert(length < len(helix))
    numberOfHelices = len(helix) - length
    for i in range(numberOfHelices + 1):
        newHelixId = '%s[%d:%d]' % (helix.id, i, length + i)
        yield AAReadWithX(newHelixId, helix.sequence[i:length + i])


allHelices = [helix for helix in FastaReads(sys.stdin,
                                            readClass=AAReadWithX,
                                            checkAlphabet=0)]

allLengths = [len(helix) for helix in allHelices]
minLength = min(allLengths)
maxLength = max(allLengths)

for length in range(minLength, maxLength + 1):
    for helix in allHelices:
        if length == len(helix):
            newHelix = AAReadWithX(helix.id, helix.sequence)
            print(newHelix.toString(format_='fasta'), file=sys.stdout)
        elif length < len(helix):
            for newHelix in makeSubstring(length, helix):
                print(newHelix.toString(format_='fasta'), file=sys.stdout)
