#!/usr/bin/env python

import sys
import argparse
from time import time
from os.path import basename

from dark.fasta import FastaReads
from dark.reads import AARead

from light.landmarks import findLandmark, ALL_LANDMARK_FINDER_CLASSES
from light.trig import findTrigPoint, ALL_TRIG_FINDER_CLASSES
from light.database import Database

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
        choices=sorted(klass.NAME for klass in ALL_LANDMARK_FINDER_CLASSES),
        help='The name of a landmark finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--trig', action='append', type=str, dest='trigFinderNames',
        choices=sorted(klass.NAME for klass in ALL_TRIG_FINDER_CLASSES),
        help='The name of a trig point finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--limitPerLandmark', type=int, default=None,
        help='A limit on the number of pairs to yield per landmark per read.')

    parser.add_argument(
        '--maxDistance', type=int, default=None,
        help='The maximum distance permitted between yielded pairs.')

    parser.add_argument(
        '--minDistance', type=int, default=None,
        help='The minimum distance permitted between yielded pairs.')

    parser.add_argument(
        '--bucketFactor', type=int, default=1,
        help=('A factor by which the distance between landmark and trig point '
              'is divided.'))

    args = parser.parse_args()

    landmarkFinderNames = (args.landmarkFinderNames or
                           [klass.NAME for klass in
                            ALL_LANDMARK_FINDER_CLASSES])
    trigFinderNames = (args.trigFinderNames or
                       [klass.NAME for klass in ALL_TRIG_FINDER_CLASSES])

    if len(landmarkFinderNames) + len(trigFinderNames) == 0:
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
            print '%s: Could not find trig point finder %r.' % (
                basename(sys.argv[0]), trigFinderName)
            sys.exit(1)

    # Create the database, add reads to it, print statistics.
    database = Database(landmarkFinderClasses, trigFinderClasses,
                        args.limitPerLandmark, args.maxDistance,
                        args.minDistance, args.bucketFactor)
    databaseReads = FastaReads(args.databaseFasta, readClass=AARead)
    for read in databaseReads:
        database.addRead(read)

    buildDbTime = time()

    print 'Database built in %.2f seconds.' % (buildDbTime - startTime)
    print 'Parameters:', database.saveParamsAsJSON()

    lookupStartTime = time()

    matches = 0
    readMatch = 0
    lookupReads = FastaReads(args.sequenceFasta, readClass=AARead)
    for read in lookupReads:
        result = database.find(read)
        if len(result.significant) > 0:
            readMatch += 1
            for subjectIndex in result.significant:
                matches += 1
                score = result.significant[subjectIndex]['matchScore']
                print 'Query: %s, Subject %s, score: %d' % (
                    read.id, database.readInfo[subjectIndex][0], score)
        else:
            print read.id, 'Read not found'

    finishedTime = time()

    print 'Look up done in %.2f seconds' % (finishedTime - lookupStartTime)
    print ('Run performance test in %.2f seconds. %d out of %d reads match. '
           'Found %d significant matches.') % (
          (finishedTime - lookupStartTime), readMatch, len(lookupReads),
        matches)
