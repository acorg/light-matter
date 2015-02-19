#!/usr/bin/env python

import argparse

from dark.fasta import combineReads

from light.database import DatabaseSpecifier


def printResult(result, database, verbose):
    """
    Print the details of trying to match a read against the database in a
    human-readable format.

    @param result: A C{light.result.Result} instance.
    @param database: A C{light.database.Database} instance.
    @param verbose: If C{True}, print details of landmark and trig point
        matches.
    """
    significant = set(result.significant())
    scannedRead = result.scannedRead
    read = scannedRead.read
    coveredIndices = len(result.scannedRead.coveredIndices())

    print 'Read: %s' % read.id
    print '  Sequence: %s' % read.sequence
    print '  Length: %d' % len(read.sequence)
    print '  Covered indices: %d (%.2f%%)' % (
        coveredIndices, coveredIndices / float(len(read.sequence)) * 100.0)

    # Print read landmarks and trig points.
    print '  Landmark count %d, trig point count %d' % (
        len(scannedRead.landmarks), len(scannedRead.trigPoints))
    if verbose:
        for landmark in scannedRead.landmarks:
            print '    ', landmark
        for trigPoint in scannedRead.trigPoints:
            print '    ', trigPoint

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
        histogram = analysis['histogram']
        binCounts = [len(b) for b in histogram.bins]
        print '    Title: %s' % title
        print '      Score: %s' % analysis['score']
        print '      Hash count: %s' % analysis['hashCount']
        print '      Sequence: %s' % sequence
        print '      Database subject index: %d' % subjectIndex
        print '      Histogram'
        print '        Number of bins: %d' % len(histogram.bins)
        print '        Significance cutoff: %s' % (
            analysis['significanceCutoff'])
        print '        Significant bin count: %s' % (
            analysis['significantBinCount'])
        print '        Max bin count: %r' % max(binCounts)
        print '        Counts: %r' % binCounts
        print '        Max offset delta: %d' % histogram.max
        print '        Min offset delta: %d' % histogram.min


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
        help='Amino acid sequences to add to the set of reads to be looked up '
        'in the database. The sequence id will be the text up to the last '
        'space, if any, otherwise will be automatically assigned.')

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

    parser.add_argument(
        '--verbose', default=False, action='store_true',
        help='If True, print details of landmark and trig point matches.')

    databaseSpecifier = DatabaseSpecifier()
    databaseSpecifier.addArgsToParser(parser)

    args = parser.parse_args()

    database = databaseSpecifier.getDatabaseFromArgs(args)
    reads = combineReads(args.fastaFile, args.sequences)

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
            printResult(result, database, args.verbose)
        else:
            result.save()
