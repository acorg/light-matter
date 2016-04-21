#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from time import time

from light.database import DatabaseSpecifier
from light.parameters import DatabaseParameters
from light.performance.query import queryDatabase

"""
REMARKS:

The test files are:
Test1: gi|526245012|ref|YP_008320338.1| terminase large subunit
       [Paenibacillus phage phiIBB_Pl23]
Test2: Same as test 1, but split onto 70aa non-overlapping chunks.
Test3: Database: Bunyaviridae, S-segments, downloaded from ViPR, 8/10/2014
       Read: gb:AFD33438|gi:380042477|UniProtKB:H9BQP9|Organism:Amur virus
             NA33|Protein Name:nucleocapsid protein|Gene Symbol:S|Segment:S
Test4: Same as test 3, but read split onto 70aa non-overlapping chunks.
"""

T1DB = 'performance/database/t1-db.fasta'
T1READ = 'performance/read/t1-read.fasta'

T2DB = 'performance/database/t2-db.fasta'
T2READ = 'performance/read/t2-read.fasta'

T3DB = 'performance/database/t3-db.fasta'
T3READ = 'performance/read/t3-read.fasta'

T4DB = 'performance/database/t4-db.fasta'
T4READ = 'performance/read/t4-read.fasta'


class WriteMarkdownFile(object):
    """
    A class that writes results of tests to a file.

    @param outputFile: a C{str} filename where the output will be written to.
    """
    def __init__(self, outputFile, verbose):
        self.outputFile = outputFile
        self.verbose = verbose

    def open(self):
        self.openedFile = open(self.outputFile, 'w')

    def writeHeader(self, dbParams):
        self.openedFile.write('Title:\nDate:\nCategory: light-matter\nTags: '
                              'light-matter, benchmarking\nSummary: '
                              'Performance and sensitivity testing\n\n'
                              '#####Database arguments:</b> '
                              'Landmarks: %s, trig points: %s, '
                              'maxDistance: %s, minDistance: %s '
                              'limitPerLandmark: %s, distanceBase %f\n\n' %
                              (dbParams.landmarkFinderNames(),
                               dbParams.trigPointFinderNames(),
                               dbParams.maxDistance, dbParams.minDistance,
                               dbParams.limitPerLandmark,
                               dbParams.distanceBase))

    def writeTest(self, testName, testResult, time, readNr):
        """
        Collect and print the test result.

        @param testName: the C{str} name of the test.
        @param testResult: a C{dict} of test results.
        @param time: the C{float} time it took to run the test.
        """
        self.openedFile.write('####%s\n\nRun in %.2f seconds.\n%d reads out '
                              'of %d matched.\n\n' % (testName, time,
                                                      len(testResult), readNr))
        for read in testResult:
            self.openedFile.write('Read %s matched %d subjects.\n\n' % (
                                  read, len(testResult[read])))
        if self.verbose:
            for read in testResult:
                for match in testResult[read]:
                    self.openedFile.write('Query: %s, Subject: %s, Score: %d '
                                          '\n\n' % (read, match,
                                                    testResult[read][match]))

    def close(self):
        self.openedFile.close()


if __name__ == '__main__':
    startTime = time()

    parser = argparse.ArgumentParser(
        description=('A script to run a set of performance and sensitivity '
                     'tests and generates a HTML file with the results.'))

    parser.add_argument(
        '--outputFile', required=True, default='index.md',
        help='The filename where the output will be written to.')

    parser.add_argument(
        '--verbose', default=False, action='store_true',
        help=('If True, information about each matching read will be written '
              'out.'))

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)

    args = parser.parse_args()
    dbParams = DatabaseParameters.fromArgs(args)

    # start writing to file
    writer = WriteMarkdownFile(args.outputFile, args.verbose)
    writer.open()
    writer.writeHeader(dbParams)

    # run tests
    # 1) A complete sequence must match itself:
    oneStart = time()
    print('1) A complete sequence must find itself.', file=sys.stderr)
    database = databaseSpecifier.getDatabaseFromArgs(args, dbParams)
    oneResult = queryDatabase(T1DB, T1READ, database)
    oneTime = time() - oneStart
    writer.writeTest('1) A complete sequence must match itself', oneResult,
                     oneTime, 1)

    # 2) Reads made from a sequence must match itself:
    twoStart = time()
    print('2) Reads made from a sequence must match itself.', file=sys.stderr)
    database = databaseSpecifier.getDatabaseFromArgs(args, dbParams)
    twoResult = queryDatabase(T2DB, T2READ, database)
    twoTime = time() - twoStart
    writer.writeTest('2) Reads made from a sequence must match itself',
                     twoResult, twoTime, 9)

    # 3) A sequence must find related sequences:
    threeStart = time()
    print('3) A sequence must match related sequences.', file=sys.stderr)
    database = databaseSpecifier.getDatabaseFromArgs(args, dbParams)
    threeResult = queryDatabase(T3DB, T3READ, database)
    threeTime = time() - threeStart
    writer.writeTest('3) A sequence must find related sequences', threeResult,
                     threeTime, 1)

    # 4) Reads made from a sequence must match related sequences:
    fourStart = time()
    print(('4) Reads made from a sequence must match related '
           'sequences.'), file=sys.stderr)
    database = databaseSpecifier.getDatabaseFromArgs(args, dbParams)
    fourResult = queryDatabase(T4DB, T4READ, database)
    fourTime = time() - fourStart
    fourTitle = '4) Reads made from a sequence must match related sequences'
    writer.writeTest(fourTitle, fourResult, fourTime, 7)

    # close file
    writer.close()
