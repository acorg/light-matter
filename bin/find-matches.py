#!/usr/bin/env python

import argparse
from time import time
from collections import defaultdict

from dark.fasta import FastaReads

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

    # read database from file.
    database = ScannedReadDatabase().load(args.databaseFile)

    # make a database for look up:
    readsDatabase = ScannedReadDatabase(database.landmarkFinderClasses,
                                        database.trigFinderClasses,
                                        database.limitPerLandmark,
                                        database.maxDistance)
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
