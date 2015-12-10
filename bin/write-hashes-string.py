#!/usr/bin/env python

from __future__ import print_function

import argparse

from light.database import DatabaseSpecifier
from light.performance.hashes import HashesString


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=('For a set of given sequences, find all hashes and for '
                     'each sequence write out a string of 1 or 0 denoting '
                     'whether a hash is present in that sequence or not. '
                     'Only include hashes if they occur in more than at least '
                     'a specified percentage of all given sequences.'))

    parser.add_argument(
        '--fastaFile',
        help=('The name of the FASTA file with sequences.'))

    parser.add_argument(
        '--cutoff', type=float, default=0.0,
        help=('Each hash has to be present in at least this (float) fraction '
              'of all sequences to be taken into account.'))

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()

    database = databaseSpecifier.getDatabaseFromArgs(args)

    hashesString = HashesString(args.fastaFile, args.cutoff, database=database)

    print(hashesString)
