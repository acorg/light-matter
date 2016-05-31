#!/usr/bin/env python

"""
Read a FASTA file from stdin, and for each sequence it contains, write all
its subsequences to stdout, also in FASTA format.

Run with --help to see all options.
"""

import sys
import argparse
from collections import Counter

from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.reads import AAReadWithX
from dark.reads import SSAAReadWithX


def makeSubsequences(length, sequence):
    """
    Take a sequence and make subsequences of a specified length.

    @param length: An C{int} length of the desired subsequences.
    @param sequence: A {dark.Read} instance (or one of its subclasses) to make
        subsequences from.
    @return: A generator that yields subsequences of C{sequence} (of the same
        type as C{sequence}). The yielded subsequences will have a [from:to]
        suffix appended to their ids if they are shorter than the original
        sequence. The C{int} from and to values follow the Python string
        indexing convention.
    """
    seqLen = len(sequence)
    if length == seqLen:
        yield sequence
    elif length <= seqLen:
        numberOfSubsequences = len(sequence) - length + 1
        for i in range(numberOfSubsequences):
            result = sequence[i:length + i]
            result.id += '[%d:%d]' % (i, length + i)
            yield result

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Read a FASTA file from stdin, and for each sequence, '
                 'write all its subsequences to stdout.'))

parser.add_argument(
    '--minLength', type=int, default=4,
    help='The minimum length subsequence to produce.')

parser.add_argument(
    '--maxLength', type=int,
    help=('The maximum length subsequence to produce. If unspecified, '
          'all subsequences at least as long as minLength will be produced.'))

parser.add_argument(
    '--pdb', default=False, action='store_true',
    help=('If True, the input is in PDB format, in which two FASTA records '
          'are used per sequence. The first holds the AA sequence and the '
          'second has the structure. If False, the input is taken as regular '
          'FASTA.'))

parser.add_argument(
    '--noSummary', default=False, action='store_true',
    help=('If specified do NOT print (to stderr) a summary showing the number '
          'of sequences read and the number of subsequences written of each '
          'different length.'))

args = parser.parse_args()

if args.pdb:
    readClass = SSAAReadWithX
    readsClass = SSFastaReads
    format_ = 'fasta-ss'
else:
    readClass = AAReadWithX
    readsClass = FastaReads
    format_ = 'fasta'

readCount = writeCount = 0
summary = Counter()

for sequence in readsClass(sys.stdin, readClass=readClass, checkAlphabet=0):
    readCount += 1
    for length in range(args.minLength, (args.maxLength or len(sequence)) + 1):
        for subsequence in makeSubsequences(length, sequence):
            summary[length] += 1
            writeCount += 1
            print(subsequence.toString(format_=format_), end='')

if not args.noSummary:
    print('Read %d sequence%s, wrote %d subsequence%s.' %
          (readCount, '' if readCount == 1 else 's',
           writeCount, '' if writeCount == 1 else 's'), file=sys.stderr)
    if writeCount:
        print('Subsequence lengths and counts:', file=sys.stderr)
        for length in sorted(summary):
            print('  ', length, summary[length], file=sys.stderr)
