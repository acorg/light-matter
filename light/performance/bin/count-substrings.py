#!/usr/bin/env python

from __future__ import print_function

from collections import Counter
import argparse
import sys

from dark.fasta_ss import SSFastaReads
from dark.reads import SSAAReadWithX


parser = argparse.ArgumentParser(
    description=('For each sequence in the specified file, make all '
                 'substrings of the length specified on stdin and count the'
                 'number of occurences for each substring.'))

parser.add_argument(
    '--length', type=int,
    help='The length of the substrings that will be produced. Must be bigger '
         'than 2.')

parser.add_argument(
    '--pdbFile', help='A filename of the pdb file containing sequences and '
    'their structure annotation.')

parser.add_argument(
    '--excludeFile', help='A filename of a file containing strings that '
    'should not be counted.')

args = parser.parse_args()


pdbReads = SSFastaReads(args.pdbFile, readClass=SSAAReadWithX,
                        checkAlphabet=0)

exclude = set()
with open(args.excludeFile) as fp:
    for line in fp:
        string = line.split(' ')[0]
        if len(string) == args.length:
            exclude.add(string)

# We only want to make subsets that are at least of length 4.
if args.length < 3:
    print('The length of the substrings must be >= 3. It\'s currently %d.' % (
          args.length), file=sys.stderr)
    sys.exit(0)

allKMers = Counter()

for i, read in enumerate(pdbReads):
    readLength = len(read.sequence)
    if args.length <= readLength:
        numberOfStrings = readLength - args.length + 1
        for i in range(numberOfStrings):
            newString = read.sequence[i:args.length + i]
            if newString not in exclude:
                allKMers[newString] += 1

for kMer in allKMers:
    print(kMer, allKMers[kMer])
