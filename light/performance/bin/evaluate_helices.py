#!/usr/bin/env python

from __future__ import print_function
import re
import sys

from dark.fasta_ss import SSFastaReads
from dark.fasta import FastaReads
from dark.reads import AAReadWithX


def evaluateMatch(structureString, start, end):
    """
    Test if a match is correct. There are four scenarios:
    1) The alpha helix matches part of a sequence that's not an alpha helix.
        --> false positive
    2) The alpha helix matches part of a sequence that's an alpha helix. The
       alpha helix in the sequence doesn't extend to the left or right.
        --> true positive.
    3) The alpha helix matches part of a sequence that's an alpha helix. The
       alpha helix in the sequence extends to the left.
        --> false positive.
    4) The alpha helix matches part of a sequence that's an alpha helix. The
       alpha helix in the sequence extends to the right.
        --> true positive.

    @param structureString: A C{str} of a structure sequence.
    @param start: An C{int} start of the match.
    @param end: An C{int} end of the match.

    @return: C{True} if the match is a true positive and C{False} if the match
        is a false positive.
    """
    assert 0 <= start < end

    if start > 0 and structureString[start - 1] == 'H':
        return False

    for aaIndex in range(start, end):
        if structureString[aaIndex] != 'H':
            return False

    return True


def evaluateHelices(helixFile, pdbFile, outFile):
    """
    For each helix in the helixFile, evaluate its true and false positives.

    @param helixFile: A C{str} filename of a file with known alpha helix
        sequences.
    @param pdbFile: A C{str} filename of the pdb file containing sequence and
        their structure annotation.
    @param outFile: A C{str} filename of the output file.
    """
    pdbReads = [(read.sequence, read.structure) for read in
                SSFastaReads(pdbFile, checkAlphabet=0)]
    helices = FastaReads(helixFile, readClass=AAReadWithX, checkAlphabet=0)

    with open(outFile, 'w') as fp:
        for i, helix in enumerate(helices):
            truePositive = falsePositive = 0
            helixSequence = helix.sequence
            if ('X' in helixSequence or 'Z' in helixSequence or 'B' in
                    helixSequence):
                continue
            else:
                uniqueRegex = re.compile(helixSequence)
                for sequence, structure in pdbReads:
                    for match in uniqueRegex.finditer(sequence):
                        start = match.start()
                        end = match.end()
                        if evaluateMatch(structure, start, end):
                            truePositive += 1
                        else:
                            falsePositive += 1

                fp.write('%s %d %d\n' % (helix.sequence, truePositive,
                                         falsePositive))


if __name__ == '__main__':
    if len(sys.argv) != 4:
        from os.path import basename
        print('Usage: %s helixFile, pdbFile, outFile' % basename(sys.argv[0]),
              file=sys.stderr)
    else:
        helixFile, pdbFile, outFile = sys.argv[1:]
        evaluateHelices(helixFile, pdbFile, outFile)
