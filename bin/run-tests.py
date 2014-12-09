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
        description='A script for performance testing. Test if a seuqnece '
        'is found in a database.')

    parser.add_argument(
        '--databaseFasta', required=True,
        help='The name of the FASTA file containing the sequences that should '
        'be turned into a database.')

    parser.add_argument(
        '--sequenceFasta', type=str, required=True,
        help='The name of the FASTA file that contains the sequence to be '
        'looked up.')

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

    # Create the database, add reads to it, print statistics.
    database = ScannedReadDatabase(landmarkFinderClasses, trigFinderClasses,
                                   args.limitPerLandmark, args.maxDistance)
    databaseReads = FastaReads(args.databaseFasta)
    for read in databaseReads:
        database.addRead(read)

    buildDbTime = time()

    print >>sys.stderr, 'Database built in %.2f seconds.' % (
        buildDbTime - startTime)
    print >>sys.stderr, 'Parameters:'
    database.saveParamsAsJSON()

    lookupStartTime = time()

    matches = 0
    readMatch = 0
    lookupReads = FastaReads(args.sequenceFasta)
    for read in lookupReads:
        result = database.find(read)
        if len(result.significant) > 0:
            readMatch += 1
            for subjectIndex in result.significant:
                matches += 1
                score = result.significant[subjectIndex]['matchScore']
                print >>sys.stderr, 'Query: %s, Subject %s, score: %d' % (
                    read.id, database.readInfo[subjectIndex][0], score)
        else:
            print >>sys.stderr, read.id, 'Read not found'

    finishedTime = time()

    print >>sys.stderr, 'Look up done in %.2f seconds' % (
        finishedTime - lookupStartTime)
    print >>sys.stderr, ('Run performance test in %.2f seconds. %d out of '
                         '%d reads match. Found %d significant matches.') % (
        (finishedTime - lookupStartTime), readMatch, len(lookupReads), matches)
