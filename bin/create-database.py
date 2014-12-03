#!/usr/bin/env python

import sys
import argparse
from time import time
from os.path import basename

from dark.fasta import FastaReads

from light.landmarks import find as findLandmark, ALL_LANDMARK_FINDER_CLASSES
from light.trig import find as findTrigPoint, ALL_TRIG_FINDER_CLASSES
from light.database import ScannedReadDatabase


if __name__ == '__main__':
    startTime = time()

    parser = argparse.ArgumentParser(
        description='Create a light-matter database from sequences in a '
        'FASTA file.')

    parser.add_argument(
        '--fastaFile', required=True,
        help='The name of the FASTA file containing the sequences for the '
        'database.')

    parser.add_argument(
        '--landmark', action='append', type=str, dest='landmarkFinderNames',
        default=[klass.NAME for klass in ALL_LANDMARK_FINDER_CLASSES],
        choices=sorted(klass.NAME for klass in ALL_LANDMARK_FINDER_CLASSES),
        help='The name of a landmark finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--trig', action='append', type=str, dest='trigFinderNames',
        default=[klass.NAME for klass in ALL_TRIG_FINDER_CLASSES],
        choices=sorted(klass.NAME for klass in ALL_TRIG_FINDER_CLASSES),
        help='The name of a trig point finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--limitPerLandmark', type=int, default=None,
        help='A limit on the number of pairs to yield per landmark per read.')

    parser.add_argument(
        '--maxDistance', type=int, default=None,
        help='The maximum distance permitted between yielded pairs.')

    args = parser.parse_args()

    if len(args.landmarkFinderNames) + len(args.trigFinderNames) == 0:
        print >>sys.stderr, ('You must specify either landmark or trig point '
                             'finders to find.\n%s') % parser.format_usage()
        sys.exit(1)

    # Make sure all landmark finders requested exist.
    landmarkFinderClasses = []
    for landmarkFinderName in args.landmarkFinderNames:
        landmarkFinderClass = findLandmark(landmarkFinderName)
        if landmarkFinderClass:
            landmarkFinderClasses.append(landmarkFinderClass)
        else:
            print >>sys.stderr, '%s: Could not find landmark finder %r.' % (
                basename(sys.argv[0]), landmarkFinderName)
            sys.exit(1)

    # Make sure all trig point finders requested exist.
    trigFinderClasses = []
    for trigFinderName in args.trigFinderNames:
        trigFinderClass = findTrigPoint(trigFinderName)
        if trigFinderClass:
            trigFinderClasses.append(trigFinderClass)
        else:
            print >>sys.stderr, '%s: Could not find trig point finder %r.' % (
                basename(sys.argv[0]), landmarkFinderName)
            sys.exit(1)

    # Create the database, add reads to it, print statistics, save to stdout.
    database = ScannedReadDatabase(landmarkFinderClasses, trigFinderClasses,
                                   args.limitPerLandmark, args.maxDistance)
    reads = FastaReads(args.fastaFile)
    for read in reads:
        database.addRead(read)
    print >>sys.stderr, database
    print >>sys.stderr, 'Database built in %.2f seconds. Saving.' % (
        time() - startTime)
    database.save()
