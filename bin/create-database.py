#!/usr/bin/env python

import sys
import argparse
from time import time

from light.database import DatabaseSpecifier


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create a light-matter database from sequences in a '
        'FASTA file and/or given on the command line.')

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()
    startTime = time()
    database = databaseSpecifier.getDatabaseFromArgs(args)
    database.print_(fp=sys.stderr)

    print('Database built in %.2f seconds. Saving...' % (
        time() - startTime), end=' ', file=sys.stderr)
    saveStartTime = time()
    database.save()
    print('saved in %.2f seconds.' % (time() - saveStartTime), file=sys.stderr)
