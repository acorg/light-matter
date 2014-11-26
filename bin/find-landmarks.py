#!/usr/bin/env python

import sys
import argparse
from os.path import basename

from dark.reads import AARead
from dark.fasta import FastaReads

from light.landmarks import find


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find landmarks in sequences in a FASTA file.')

    parser.add_argument(
        '--fastaFile', required=True,
        help='The name of the FASTA file to read.')

    parser.add_argument(
        '--landmark', required=True,
        help='The name of the landmark to look for.')

    args = parser.parse_args()

    reads = FastaReads(args.fastaFile, AARead)
    landmarkFinderClass = find(args.landmark)

    if not landmarkFinderClass:
        print >>sys.stderr, '%s: Could not find landmark class %r.' % (
            basename(sys.argv[0]), args.landmark)
        sys.exit(1)

    landmarkFinder = landmarkFinderClass().find

    for read in reads:
        seenLandmark = False
        for landmark in landmarkFinder(read):
            if seenLandmark is False:
                print read.id
                seenLandmark = True
            print '\t', landmark
