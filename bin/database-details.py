#!/usr/bin/env python

import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from light.database import Database
from light.landmarks import ALL_LANDMARK_CLASSES
from light.trig import ALL_TRIG_CLASSES

_SYMBOL_NAME = dict((finder.SYMBOL, finder.NAME) for finder in
                    ALL_LANDMARK_CLASSES | ALL_TRIG_CLASSES)


def getFinderNameFromSymbol(symbol):
    """
    Get the name of a finder class given a hashkey that was
    created for it.

    @param symbol: A C{str} symbol for the class.
    @return: A C{str} symbol name.
    """
    # Look for the entire symbol in _SYMBOL_NAME and if that's not present,
    # shorten it one character at a time from the right. This allows us to
    # find the correct finder class for symbols like PS00342 (the Prosite
    # class, whose symbol is PS).
    while symbol:
        try:
            return _SYMBOL_NAME[symbol]
        except KeyError:
            symbol = symbol[:-1]


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

    # Print basic database information.
    print 'Database file name: %s' % args.json
    if database.landmarkFinders:
        print 'Landmark finders:'
        print '  ' + '\n  '.join(sorted(
            finder.NAME for finder in database.landmarkFinders))
    else:
        print 'Landmark finders: none'
    if database.trigPointFinders:
        print 'Trig point finders:'
        print '  ' + '\n  '.join(sorted(
            finder.NAME for finder in database.trigPointFinders))
    else:
        print 'Trig point finders: none'
    print 'Sequences: %s' % database.subjectCount
    print 'Hashes: %d' % len(database.d)
    print 'Residues: %d' % database.totalResidues
    print 'Coverage: %.2f%%' % (float(database.totalCoveredResidues) /
                                database.totalResidues * 100.0)
    print 'Checksum: %s' % database.checksum
    print 'Number of landmark and trigpoints found:'

    if database.d:
        # Print hash counts.
        landmarksTrigpoints = defaultdict(int)
        for key in database.d:
            lm, tp, dist = key.split(':')
            landmarksTrigpoints[lm] += 1
            landmarksTrigpoints[tp] += 1

        for f in landmarksTrigpoints:
            print '  %s (%s): %d' % (getFinderNameFromSymbol(f), f,
                                     landmarksTrigpoints[f])

    # Write subjects to file.
    if args.printSubjects:
        subjects = []
        for sb in database.subjectInfo:
            subjects.append(SeqRecord(Seq(sb[1]), id=sb[0], description=''))
        with open(args.printSubjects, 'w') as fp:
            SeqIO.write(subjects, fp, 'fasta')
        print 'Wrote %d subjects to %s' % (len(database.subjectInfo),
                                           args.printSubjects)
    # Print hashes.
    if args.printHashes:
        for key, subjects in database.d.iteritems():
            print key
            for subject in subjects:
                index = subject['subjectIndex']
                print '    %s' % database.subjectInfo[index][0]
