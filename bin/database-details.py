#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
        subjects = []
        for sb in database.subjectInfo:
            subjects.append(SeqRecord(Seq(sb[1]), id=sb[0], description=''))
        with open(args.printSubjects, 'w') as fp:
            SeqIO.write(subjects, fp, 'fasta')
        print 'Wrote %d subjects to %s' % (len(database.subjectInfo),
                                           args.printSubjects)
