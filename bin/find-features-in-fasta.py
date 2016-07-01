#!/usr/bin/env python

"""
Read FASTA from stdin and examine it for features (landmarks and trig points).
Write a summary of the number of features that are found in the input for each
of the feature finders used.

The idea here is that you can provide input that contains known features (e.g.,
taken from PDB secondary structures) and check to see if our finders can find
the features.
"""

from __future__ import print_function

import sys
import argparse
from collections import Counter

from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.reads import AARead, SSAARead

from light.database import DatabaseParameters


parser = argparse.ArgumentParser(description='Find features in FASTA.')

parser.add_argument(
    '--margin', type=int, default=0,
    help=('The number of AA residues that were kept on either side of the '
          'features in the input. When checking to see if a feature was '
          'found, features that are contained entirely in either margin are '
          'not considered as matches.'))

parser.add_argument(
    '--minSequenceLength', '--msl', type=int, default=0,
    help=('Sequences in the input whose length is less than this value '
          'will not be examined for features. Note that this is value '
          'should include the margin width, if any.'))

parser.add_argument(
    '--pdb', default=False, action='store_true',
    help=('If True, the input is in PDB format, in which two FASTA records '
          'are used per sequence. The first holds the AA sequence and the '
          'second has the structure. If False, the input is taken as regular '
          'FASTA.'))

parser.add_argument(
    '--verbose', default=False, action='store_true',
    help=('If True, print detailed information about input sequences and '
          'matching features.'))

DatabaseParameters.addArgsToParser(parser)
args = parser.parse_args()
dbParams = DatabaseParameters.fromArgs(args)

landmarkFinders = dbParams.landmarkFinders
trigPointFinders = dbParams.trigPointFinders
verbose = args.verbose

if verbose:
    if landmarkFinders:
        print('Using landmark finder%s: %s' % (
            '' if len(landmarkFinders) == 1 else 's',
            ', '.join(finder.NAME for finder in landmarkFinders)))

    if trigPointFinders:
        print('Using trig point finder%s: %s' % (
            '' if len(trigPointFinders) == 1 else 's',
            ', '.join(finder.NAME for finder in trigPointFinders)))

# Note that inexactMatches might count multiple features found in the same
# read (even though the read is supposedly just one alpha helix feature,
# our finders may not see it that way - assuming they see anything at all).

exactMatches = Counter()
inexactMatches = Counter()
readCount = 0

margin = args.margin
minSequenceLength = args.minSequenceLength

if args.pdb:
    readClass = SSAARead
    readsClass = SSFastaReads
else:
    readClass = AARead
    readsClass = FastaReads

for read in readsClass(sys.stdin, readClass=readClass, checkAlphabet=0):
    readLen = len(read)

    if readLen < minSequenceLength:
        continue

    readCount += 1
    featureLengthInRead = readLen - 2 * margin

    if verbose:
        if readClass is AARead:
            readInfo = 'id=%s seq=%r' % (read.id, read.sequence)
        else:
            readInfo = 'id=%s seq=%r struct=%r' % (
                read.id, read.sequence, read.structure)
        print('read %d: featureLen=%d %s' % (
            readCount, featureLengthInRead, readInfo))

    rightMarginOffset = readLen - margin
    for finder in landmarkFinders + trigPointFinders:
        finderName = finder.NAME
        for feature in finder.find(read):
            if (feature.offset + feature.length < margin or
                    feature.offset >= rightMarginOffset):
                # A feature was found, but it's entirely in the left margin
                # or right margin of the read, so we ignore it.
                match = 'margin'
            elif (feature.offset == margin and
                  feature.length == featureLengthInRead):
                exactMatches[finderName] += 1
                match = 'exact'
            else:
                inexactMatches[finderName] += 1
                match = 'inexact'
            if verbose:
                print('    %s (%s) at %d length %d' % (
                    feature.name, match, feature.offset, feature.length))

print('Processed %d reads' % readCount)

for finder in landmarkFinders + trigPointFinders:
    finderName = finder.NAME
    print('%s: %d exact matches (%.2f%%), %d inexact matches (%.2f%%)' % (
        finderName, exactMatches[finderName],
        exactMatches[finderName] / readCount * 100.0,
        inexactMatches[finderName],
        inexactMatches[finderName] / readCount * 100.0))
