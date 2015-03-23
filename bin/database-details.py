#!/usr/bin/env python

import argparse

from light.database import DatabaseSpecifier


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Print information about a database.')

    parser.add_argument(
        '--printSubjects', metavar='SUBJECT-FILENAME', default=False, type=str,
        help='If a filename is given, write all subjects to the file in fasta '
        'format.')

    parser.add_argument(
        '--printHashes', default=False, action='store_true',
        help='If True, print all hashes and associated subjects.')

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)
    args = parser.parse_args()
    database = databaseSpecifier.getDatabaseFromArgs(args)

    database.print_(printHashes=args.printHashes)

    # Write subjects to file.
    if args.printSubjects:
        with open(args.printSubjects, 'w') as fp:
            for subject in database.getSubjects():
                print >>fp, subject.toString('fasta')
        print 'Wrote %d subjects to %s' % (len(database.subjectCount),
                                           args.printSubjects)
