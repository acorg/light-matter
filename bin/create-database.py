#!/usr/bin/env python

import sys
import argparse
from time import time
from os.path import basename

from dark.fasta import FastaReads
from dark.reads import Reads, AARead

from light.landmarks import find as findLandmark, ALL_LANDMARK_FINDER_CLASSES
from light.trig import find as findTrigPoint, ALL_TRIG_FINDER_CLASSES
from light.database import Database


if __name__ == '__main__':
    startTime = time()

    parser = argparse.ArgumentParser(
        description='Create a light-matter database from sequences in a '
        'FASTA file.')

    parser.add_argument(
        '--fastaFile',
        help='The name of the FASTA file containing the sequences for the '
        'database.')

    parser.add_argument(
        '--sequence', action='append', dest='sequences',
        metavar='"id sequence"',
        help='Amino acid sequences to add to the database. The sequence id '
        'will be the text up to the last space, if any, otherwise will be '
        'automatically assigned.')

    parser.add_argument(
        '--landmark', action='append', dest='landmarkFinderNames',
        choices=sorted(klass.NAME for klass in ALL_LANDMARK_FINDER_CLASSES),
        help='The name of a landmark finder to use. May be specified '
        'multiple times.')

    parser.add_argument(
        '--trig', action='append', dest='trigFinderNames',
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
    for landmarkFinderName in landmarkFinderNames:
        landmarkFinderClass = findLandmark(landmarkFinderName)
        if landmarkFinderClass:
            landmarkFinderClasses.append(landmarkFinderClass)
        else:
            print >>sys.stderr, '%s: Could not find landmark finder %r.' % (
                basename(sys.argv[0]), landmarkFinderName)
            sys.exit(1)

    # Make sure all trig point finders requested exist.
    trigFinderClasses = []
    for trigFinderName in trigFinderNames:
        trigFinderClass = findTrigPoint(trigFinderName)
        if trigFinderClass:
            trigFinderClasses.append(trigFinderClass)
        else:
            print >>sys.stderr, '%s: Could not find trig point finder %r.' % (
                basename(sys.argv[0]), trigFinderName)
            sys.exit(1)

    # Create the database, add reads to it, print statistics, save to stdout.
    database = Database(landmarkFinderClasses, trigFinderClasses,
                        args.limitPerLandmark, args.maxDistance)

    # Arrange to read AA subject sequences from a FASTA file, if given.
    if args.fastaFile:
        reads = FastaReads(args.fastaFile)
    else:
        if args.sequences is None:
            print >>sys.stderr, (
                'No reads to search for given.\n%s' % parser.format_usage()),
            sys.exit(1)
        reads = Reads()

    # Plus any manually-specified AA subject sequences.
    if args.sequences:
        for count, sequence in enumerate(args.sequences, start=1):
            # Try splitting the sequence on its last space and using the first
            # part of the split as the read id. If there is no space, assign a
            # generic id.  Note that using 'command-line-read-' as the read id
            # prefix could collide with ids in the FASTA file, if given. So
            # output might be ambiguous. But that seems pretty unlikely.
            parts = sequence.rsplit(' ', 1)
            if len(parts) == 2:
                readId, sequence = parts
            else:
                readId = 'command-line-read-%d' % count
            reads.add(AARead(readId, sequence))

    for read in reads:
        database.addSubject(read)

    print >>sys.stderr, database
    print >>sys.stderr, 'Database built in %.2f seconds. Saving...' % (
        time() - startTime),

    saveStartTime = time()
    database.save()
    print >>sys.stderr, 'saved in %.2f seconds.' % (time() - saveStartTime)
