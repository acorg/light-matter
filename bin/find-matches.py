#!/usr/bin/env python

import argparse

from dark.fasta import combineReads

from light.database import DatabaseSpecifier


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
            result.print_(database, verbose=args.verbose)
        else:
            result.save()
