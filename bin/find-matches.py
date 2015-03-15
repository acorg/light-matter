#!/usr/bin/env python

import sys
from os.path import basename
import argparse

from dark.fasta import combineReads

from light.database import DatabaseSpecifier


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Look up reads in a database.')

    parser.add_argument(
        '--fastaFile',
        help=('The name of the FASTA file containing the sequences that '
              'should be looked up.'))

    parser.add_argument(
        '--sequence', action='append', dest='sequences',
        metavar='"id sequence"',
        help=('Amino acid sequences to add to the set of reads to be looked '
              'up in the database. The sequence id will be the text up to the '
              'first space, if any, otherwise will be automatically '
              'assigned.'))

    parser.add_argument(
        '--significanceFraction', type=float, default=None,
        help=('The (float) fraction of all (landmark, trig point) pairs for a '
              'scannedRead that need to fall into the same histogram bucket '
              'for that bucket to be considered a significant match with a '
              'database title.'))

    parser.add_argument(
        '--human', default=False, action='store_true',
        help=('If True, produce human-readable output (this is useful when '
              'informally examining sequences from the command line). Do not '
              'print database parameters, and show detailed information about '
              'matches.'))

    parser.add_argument(
        '--printHistograms', default=False, action='store_true',
        help=('If True, and --human is specified, print details of histogram '
              'bin counts.'))

    parser.add_argument(
        '--printSequences', default=False, action='store_true',
        help=('If True, and --human is specified, print query and subject '
              'sequences.'))

    parser.add_argument(
        '--excludeFeatures', default=False, action='store_true',
        help=('If True, and --human is specified, print details of histogram '
              'bin counts.'))

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)

    args = parser.parse_args()

    # Check for argument compatibility. It would be more friendly to
    # automatically turn on --human if any of --printHistograms or
    # --printSequences or --excludeFeatures is used, but I think that could
    # be confusing, so for now let's be strict.
    if ((args.printHistograms or args.printSequences or args.excludeFeatures)
            and not args.human):
        print >>sys.stderr, (
            '%s: If you specify --printHistograms or --printSequences or '
            '--excludeFeatures, you must also use --human.'
            % basename(sys.argv[0]))
        sys.exit(1)

    database = databaseSpecifier.getDatabaseFromArgs(args)
    reads = combineReads(args.fastaFile, args.sequences)
    significanceFraction = args.significanceFraction

    # Look up each read in the database and print its matches, either in
    # human-readable form or as JSON.
    if args.human:
        for count, read in enumerate(reads, start=1):
            result = database.find(
                read, significanceFraction=significanceFraction,
                storeFullAnalysis=True)

            result.print_(database,
                          printSequences=args.printSequences,
                          printFeatures=not args.excludeFeatures,
                          printHistograms=args.printHistograms,
                          queryDescription='Query %d title' % count)
            if count > 1:
                print '---'
    else:
        database.saveParamsAsJSON()
        for read in reads:
            result = database.find(
                read, significanceFraction=significanceFraction)
            result.save()
