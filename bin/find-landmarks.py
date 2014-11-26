#!/usr/bin/env python

import sys
import argparse
from os.path import basename

from dark.reads import AARead
from dark.fasta import FastaReads

from light.landmarks import find, ALL_LANDMARK_FINDER_CLASSES


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find landmarks in sequences in a FASTA file.')

    parser.add_argument(
        '--fastaFile', required=True,
        help='The name of the FASTA file to read.')

    parser.add_argument(
        '--landmark', required=True, action='append', type=str,
        choices=sorted(c.__name__ for c in ALL_LANDMARK_FINDER_CLASSES),
        help='The name of the landmark finder to use. May be specified '
        'multiple times.')

    args = parser.parse_args()

    landmarkFinders = []
    for landmarkFinderName in args.landmark:
        landmarkFinderClass = find(landmarkFinderName)
        if landmarkFinderClass:
            landmarkFinders.append(landmarkFinderClass().find)
        else:
            print >>sys.stderr, '%s: Could not find landmark finder %r.' % (
                basename(sys.argv[0]), landmarkFinderName)
            sys.exit(1)

    # For each read, find all landmarks.
    for read in FastaReads(args.fastaFile, AARead):
        seenAnyLandmark = False
        for landmarkFinder in landmarkFinders:
            for landmark in landmarkFinder(read):
                if not seenAnyLandmark:
                    print read.id
                    seenAnyLandmark = True
                print '\t', landmark
