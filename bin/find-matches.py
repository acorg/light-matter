#!/usr/bin/env python

import sys
import argparse
from time import time
from os.path import basename

from dark.fasta import FastaReads

from light.landmarks import find as findLandmark, ALL_LANDMARK_FINDER_CLASSES
from light.trig import find as findTrigPoint, ALL_TRIG_FINDER_CLASSES
from light.database import ScannedReadDatabase, evaluate

if __name__ == '__main__':
    startTime = time()

    parser = argparse.ArgumentParser(
        description='Look up reads in a database.')

    parser.add_argument(
        '--fastaFile', required=True,
        help='The name of the FASTA file containing the sequences that should '
        'be looked up.')

    parser.add_argument(
        '--databaseFile', type=str, required=True,
        help='The name of the file that contains the database.')

    args = parser.parse_args()

    # read out trigpoints, landmarks, maxDistance, limitperlandmark out of
    # database file.
    #database = ScannedReadDatabase().read(args.databaseFile)
    #--> function that terry is writing now.
    #landmarkFinderNames =
    #trigFinderNames =
    #limitPerLandmark =
    #maxDistance =
    #database=

    # make a database for look up:
    readsDatabase = ScannedReadDatabase(landmarkFinderClasses,
                                        trigFinderClasses, limitPerLandmark,
                                        maxDistance)
    reads = FastaReads(args.fastaFile)
    found = defaultdict(defaultdict(list))
    for read in reads:
        for result in readsDatabase.find(read):
            found[result.subjectId][result.queryId].append(result.offset)

    # evaluate the significance of the found matches.
    significant = evaluate(found)

    print 'Significant matches found:'
    for match in significant:
        print 'Subject: %s, Query: %s, count: %d' % (match[0], match[1],
                                                     match[2])
