#!/usr/bin/env python

import argparse

from dark.fasta import FastaReads

from light.database import ScannedReadDatabase

if __name__ == '__main__':

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

    # Read database from file.
    with open(args.databaseFile) as fp:
        database = ScannedReadDatabase.load(fp)

    database.saveParamsAsJSON()

    reads = FastaReads(args.fastaFile)
    for read in reads:
        result = database.find(read)
        result.save()
