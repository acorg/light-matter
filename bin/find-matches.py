#!/usr/bin/env python

import sys
import argparse

from dark.fasta import FastaReads
from dark.reads import Reads, AARead

from light.database import Database


def printResult(result, database):
    """
    Print the details of trying to match a read against the database in a
    human-readable format.

    @param result: A C{light.result.Result} instance.
    @param database: A C{light.database.Database} instance.
    """
    significant = set(result.significant())

    print 'Read: %s' % read.id
    print '  Sequence: %s' % read.sequence
    print '  Length: %d' % len(read.sequence)
    print '  Significant matches: %d' % len(significant)
    print '  Overall matches: %d' % len(result.matches)
    if result.matches:
        print '  Subject matches:'

    # Print matched subjects in order of decreasing score.
    subjectIndices = sorted(result.matches.keys(), reverse=True,
                            key=lambda index: result.analysis[index]['score'])
    for subjectIndex in subjectIndices:
        title, sequence = database.subjectInfo[subjectIndex]
        analysis = result.analysis[subjectIndex]
        print '    Title: %s' % title
        print '      Score: %s' % analysis['score']
        print '      Significant: %s' % analysis['significant']
        print '      Sequence: %s' % sequence
        print '      Offsets: %s' % analysis['offsets']
        print '      Histogram: %s' % analysis['histogram']
        print '      Histogram buckets: %s' % analysis['histogramBuckets']
        print '      Max bucket count: %d' % analysis['maxCount']
        print '      Mean bucket count: %.3f' % analysis['meanCount']
        print '      Database subject index: %d' % subjectIndex


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Look up reads in a database.')

    parser.add_argument(
        '--fastaFile',
        help='The name of the FASTA file containing the sequences that should '
        'be looked up.')

    parser.add_argument(
        '--sequence', action='append', dest='sequences',
        metavar='"id sequence"',
        help='Amino acid sequences to add to the database. The sequence id '
        'will be the text up to the last space, if any, otherwise will be '
        'automatically assigned.')

    parser.add_argument(
        '--databaseFile', required=True,
        help='The name of the file that contains the database.')

    parser.add_argument(
        '--significanceFraction', type=float, default=None,
        help='The (float) fraction of all (landmark, trig point) pairs for a '
        'scannedRead that need to fall into the same histogram bucket for '
        'that bucket to be considered a significant match with a database '
        'title.')

    parser.add_argument(
        '--human', default=False, action='store_true',
        help='If True, produce human-readable output (this is useful when '
        'informally examining sequences from the command line). Do not print '
        'database parameters, and show detailed information about matches.')

    args = parser.parse_args()

    # Add AA sequences from a FASTA file, if given.
    if args.fastaFile:
        reads = FastaReads(args.fastaFile, readClass=AARead)
    else:
        if args.sequences is None:
            print >>sys.stderr, (
                'No reads to search for given.\n%s' % parser.format_usage()),
            sys.exit(1)
        reads = Reads()

    # Add any manually-specified AA read sequences.
    if args.sequences:
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

    # Read database from file.
    with open(args.databaseFile) as fp:
        database = Database.load(fp)

    human = args.human
    if not human:
        database.saveParamsAsJSON()

    # Look up each read in the database.
    for count, read in enumerate(reads):
        result = database.find(read, args.significanceFraction,
                               storeFullAnalysis=human)
        if human:
            if count:
                print '---'
            printResult(result, database)
        else:
            result.save()
