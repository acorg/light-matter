#!/usr/bin/env python

import argparse
from operator import attrgetter
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from light.database import Database
from light.landmarks import ALL_LANDMARK_FINDER_CLASSES
from light.trig import ALL_TRIG_FINDER_CLASSES


symbolName = {}
key = attrgetter('SYMBOL')
for finder in sorted(ALL_LANDMARK_FINDER_CLASSES, key=key):
    symbolName[finder.SYMBOL] = finder.NAME
for finder in sorted(ALL_TRIG_FINDER_CLASSES, key=key):
    symbolName[finder.SYMBOL] = finder.NAME


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Print information about a database.')

    parser.add_argument(
        'json', metavar='DATABASE-JSON-file', type=str,
        help='the JSON file of a light matter database.')

    parser.add_argument(
        '--printSubjects', metavar='SUBJECT-FILENAME', default=False, type=str,
        help='If a filename is given, write all subjects to the file in fasta '
        'format.')

    parser.add_argument(
        '--printHashes', default=False, action='store_true',
        help='If True, print all hashes and associated subjects.')

    args = parser.parse_args()

    with open(args.json) as fp:
        database = Database.load(fp)

    # print basic database information
    print 'Database file name: %s' % args.json
    print 'Sequences: %s' % database.subjectCount
    print 'Hashes: %d' % len(database.d)
    print 'Residues: %d' % database.totalResidues
    print 'Coverage: %.2f%%' % (float(database.totalCoveredResidues) /
                                database.totalResidues * 100.0)
    print 'Checksum: %s' % database.checksum()
    print 'Number of landmark and trigpoints found:'

    # print hashes and subjects
    landmarksTrigpoints = defaultdict(int)
    for key in database.d:
        lm, tp, dist = key.split(':')
        landmarksTrigpoints[lm] += 1
        landmarksTrigpoints[tp] += 1

    for f in landmarksTrigpoints:
        print '%s (%s): %d' % (symbolName[f[0]], f, landmarksTrigpoints[f])

    # write subjects to file:
    if args.printSubjects:
        subjects = []
        for sb in database.subjectInfo:
            subjects.append(SeqRecord(Seq(sb[1]), id=sb[0], description=''))
        with open(args.printSubjects, 'w') as fp:
            SeqIO.write(subjects, fp, 'fasta')
        print 'Written %d subjects to %s' % (len(database.subjectInfo),
                                             args.printSubjects)
    # print hashes
    if args.printHashes:
        for key, subjects in database.d.iteritems():
            print key
            for subject in subjects:
                index = subject['subjectIndex']
                print '    %s' % database.subjectInfo[index][0]
