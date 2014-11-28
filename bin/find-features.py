#!/usr/bin/env python

import sys
import argparse
from os.path import basename

from dark.fasta import FastaReads

from light.landmarks import find as findLandmark, ALL_LANDMARK_FINDER_CLASSES
from light.trig import find as findTrigPoint, ALL_TRIG_FINDER_CLASSES
from light.reads import ScannedRead


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find landmarks in sequences in a FASTA file.')

    parser.add_argument(
        '--fastaFile', required=True,
        help='The name of the FASTA file to read.')

    parser.add_argument(
        '--landmark', action='append', type=str,
        choices=sorted(c.__name__ for c in ALL_LANDMARK_FINDER_CLASSES),
        help='The name of the landmark finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--trig', action='append', type=str,
        choices=sorted(c.__name__ for c in ALL_TRIG_FINDER_CLASSES),
        help='The name of the trig point finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--verbose', default=False, action='store_true',
        help='If True, print details of landmark and trig point matches.')

    args = parser.parse_args()

    landmarks = args.landmark or []
    trigs = args.trig or []

    if len(landmarks) + len(trigs) == 0:
        print >>sys.stderr, ('You must specify either landmarks or trig '
                             'points to find.\n%s') % args.format_usage()
        sys.exit(1)

    # Make sure all landmark finders requested exist.
    landmarkFinders = []
    for landmarkFinderName in landmarks:
        landmarkFinderClass = findLandmark(landmarkFinderName)
        if landmarkFinderClass:
            landmarkFinders.append(landmarkFinderClass().find)
        else:
            print >>sys.stderr, '%s: Could not find landmark finder %r.' % (
                basename(sys.argv[0]), landmarkFinderName)
            sys.exit(1)

    # Make sure all trig point finders requested exist.
    trigFinders = []
    for trigFinderName in trigs:
        trigFinderClass = findTrigPoint(trigFinderName)
        if trigFinderClass:
            trigFinders.append(trigFinderClass().find)
        else:
            print >>sys.stderr, '%s: Could not find trig point finder %r.' % (
                basename(sys.argv[0]), landmarkFinderName)
            sys.exit(1)

    # For each read, find all landmarks and trig points and print details
    # of what we found.
    for read in FastaReads(args.fastaFile):
        scannedRead = ScannedRead(read)

        for landmarkFinder in landmarkFinders:
            for landmark in landmarkFinder(read):
                scannedRead.landmarks.append(landmark)

        for trigFinder in trigFinders:
            for trigPoint in trigFinder(read):
                scannedRead.trigPoints.append(trigPoint)

        if scannedRead.landmarks or scannedRead.trigPoints:
            print 'Read %s %d landmarks, %d trig points' % (
                read.id, len(scannedRead.landmarks),
                len(scannedRead.trigPoints))
            coveredIndices = len(scannedRead.coveredIndices())
            readLen = len(read)
            print '\tCoverage %d residues of %d (%.2f%%)' % (
                coveredIndices, readLen,
                float(coveredIndices) / readLen * 100.0)
            if args.verbose:
                for landmark in scannedRead.landmarks:
                    print '\t', landmark
                for trigPoint in scannedRead.trigPoints:
                    print '\t', trigPoint
