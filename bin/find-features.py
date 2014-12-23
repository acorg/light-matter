#!/usr/bin/env python

import sys
import argparse
from os.path import basename

from dark.fasta import FastaReads
from dark.reads import Reads, AARead

from light.landmarks import find as findLandmark, ALL_LANDMARK_FINDER_CLASSES
from light.trig import find as findTrigPoint, ALL_TRIG_FINDER_CLASSES
from light.reads import ScannedRead


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find landmarks in sequences in a FASTA file.')

    parser.add_argument(
        '--fastaFile',
        help='The name of the FASTA file to read and whose amino acid '
        'sequences will be searched for in the database.')

    parser.add_argument(
        '--sequence', action='append', dest='sequences',
        metavar='"id sequence"',
        help='Amino acid sequences to add to the database. The sequence id '
        'will be the text up to the last space, if any, otherwise will be '
        'automatically assigned.')

    parser.add_argument(
        '--landmark', action='append', dest='landmarks',
        choices=sorted(klass.NAME for klass in ALL_LANDMARK_FINDER_CLASSES),
        help='The name of a landmark finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--trig', action='append', dest='trigs',
        choices=sorted(klass.NAME for klass in ALL_TRIG_FINDER_CLASSES),
        help='The name of a trig point finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--verbose', default=False, action='store_true',
        help='If True, print details of landmark and trig point matches.')

    args = parser.parse_args()

    landmarks = args.landmarks or []
    trigs = args.trigs or []

    if len(landmarks) + len(trigs) == 0:
        print >>sys.stderr, ('You must specify either landmarks or trig '
                             'points to find.\n%s') % parser.format_usage(),
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
                basename(sys.argv[0]), trigFinderName)
            sys.exit(1)

    # Add AA sequences from a FASTA file, if given.
    if args.fastaFile:
        reads = FastaReads(args.fastaFile)
    else:
        if args.sequences is None:
            print >>sys.stderr, (
                'No sequences given.\n%s' % parser.format_usage()),
            sys.exit(1)
        reads = Reads()

    if args.sequences:
        # Add any manually-specified AA read sequences.
        for count, sequence in enumerate(args.sequences, start=1):
            # Try splitting the sequence on its last space and using the first
            # part of the split as the read id. If there is no space, assign a
            # generic id.  Note that using 'command-line-read-' as the read id
            # prefix could collide with ids in the FASTA file, if given. So
            # output might be ambiguous. But that seems pretty unlikely.
            parts = sequence.rsplit(' ', 1)
            if len(parts) > 1:
                readId, sequence = parts
            else:
                readId = 'command-line-read-%d' % count
            reads.add(AARead(readId, sequence))

    # For each read, find all landmarks and trig points and print details
    # of what we found.
    for read in reads:
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
        else:
            print 'Read %s, no landmarks or trig points found.' % read.id
