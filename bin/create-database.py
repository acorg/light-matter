#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from time import time

from light.database import DatabaseSpecifier
from light.parameters import DatabaseParameters


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Create a light-matter database from sequences in a '
                     'FASTA file and/or given on the command line.'))

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()
    if args.filePrefix is None:
        parser.print_help()
        raise RuntimeError('You must supply a database save file prefix.')

    startTime = time()
    dbParams = DatabaseParameters.fromArgs(args)
    database = databaseSpecifier.getDatabaseFromArgs(args, dbParams)
    print(database.print_(), file=sys.stderr)
    print('Database built in %.2f seconds. Saving...' % (
        time() - startTime), end=' ', file=sys.stderr)

    # Save.
    saveStartTime = time()
    database.save()
    print('saved in %.2f seconds.' % (time() - saveStartTime), file=sys.stderr)
