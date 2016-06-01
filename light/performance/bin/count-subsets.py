#!/usr/bin/env python

from __future__ import print_function

from collections import defaultdict
import argparse
import sys

from dark.fasta_ss import SSFastaReads
from dark.reads import SSAAReadWithX


parser = argparse.ArgumentParser(
    description=('For each helix in the given on stdin, evaluate its true and '
                 'false positives and print the result to stdout.'))

parser.add_argument(
    '--pdbFile', help='A filename of the pdb file containing sequences and '
    'their structure annotation.')

args = parser.parse_args()


pdbReads = [(read.sequence, read.structure) for read in
            SSFastaReads(args.pdbFile, readClass=SSAAReadWithX,
                         checkAlphabet=0)]

length = int(sys.stdin)

# we only want to make subsets that are at least of length 4.
if length < 4:
    sys.exit(0)

allKMers = defaultdict(int)

for i, read in enumerate(pdbReads):
    readLength = len(read.sequence)
    if length <= readLength:
        numberOfStrings = readLength - length + 1
        for i in range(numberOfStrings):
            allKMers[read.sequence[i:length + i]] += 1

for kMer in allKMers:
    print(kMer, allKMers[kMer])
