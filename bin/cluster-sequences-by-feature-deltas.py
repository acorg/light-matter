#!/usr/bin/env python

import argparse

from dark.fasta import combineReads

from light.database import DatabaseSpecifier
from light.performance.feature import FeatureEvaluator


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=('Cluster sequences based on offset distances between '
                     'their features.'))

    parser.add_argument(
        '--fastaFile',
        help=('The name of the FASTA file containing the sequences that '
              'should be clustered.'))

    parser.add_argument(
        '--labelFile',
        help=('The name of a file containing labels for reads. Each line must '
              'have a label, a space, then a read id.'))

    parser.add_argument(
        '--sequence', action='append', dest='sequences',
        metavar='"id sequence"',
        help=('Amino acid sequences to add to the set of reads to be '
              'clustered. The sequence id will be the text up to the last '
              'space, if any, otherwise will be automatically assigned.'))

    parser.add_argument(
        '--label', action='append', dest='labels',
        metavar='"label read-id"',
        help=('The label for a read. Format is label, then a single space '
              'then the read id (which must match the read id in the fasta '
              'file or a read given on the command line via --sequence.'))

    parser.add_argument(
        '--defaultLabel', default=0,
        help=('The default label to use for read ids that are not explicitly '
              'labeled.'))

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()
    database = databaseSpecifier.getDatabaseFromArgs(args)
    reads = combineReads(args.fastaFile, args.sequences)
    labels = {}

    # Process labels from --label arguments.
    if args.labels:
        for labelArg in args.labels:
            fields = labelArg.split(None, 1)
            if len(fields) != 2:
                raise ValueError('Bad label %r. --label arguments must have '
                                 'a label, a space, then a read id.' %
                                 labelArg)
            labels[fields[1]] = int(fields[0])

    # Read labels from a file, if --labelFile is given.
    if args.labelFile:
        with open(args.labelFile) as fp:
            for count, line in enumerate(fp, start=1):
                line = line[:-1]
                fields = line.split(None, 1)
                if len(fields) != 2:
                    raise ValueError('Bad label line %r. Label file %r, line '
                                     '%d. Label lines must have a label, a '
                                     'space, then a read id.' %
                                     (line, args.labelFile, count))
                labels[fields[1]] = int(fields[0])

    featureEvaluator = FeatureEvaluator(reads, labels, args.defaultLabel,
                                        database=database)

    featureEvaluator.print_()
