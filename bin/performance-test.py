#!/usr/bin/env python

import sys
import argparse
from time import time
from scipy import stats
import matplotlib.pyplot as plt
from collections import defaultdict

from light.database import DatabaseSpecifier
from light.performance.query import queryDatabase

"""
REMARKS:
Z-scores are from Cerny et al., 2014, 'Evolution of Tertiary Structure of
Viral RNA Dependent Polymerases.'
The test files are:
Test1: gi|526245012|ref|YP_008320338.1| terminase large subunit
       [Paenibacillus phage phiIBB_Pl23]
Test2: Same as test 1, but split onto 70aa non-overlapping chunks.
Test3: Database: Bunyaviridae, S-segments, downloaded from ViPR, 8/10/2014
       Read: gb:AFD33438|gi:380042477|UniProtKB:H9BQP9|Organism:Amur virus
             NA33|Protein Name:nucleocapsid protein|Gene Symbol:S|Segment:S
Test4: Same as test 3, but read split onto 70aa non-overlapping chunks.
Test5 and Test6: The polymerase sequences described in Cerny et al., 2014.
"""

# GLOBALS:
# filenames
T1DB = 'performance/database/t1-db.fasta'
T1READ = 'performance/read/t1-read.fasta'

T2DB = 'performance/database/t2-db.fasta'
T2READ = 'performance/read/t2-read.fasta'

T3DB = 'performance/database/t3-db.fasta'
T3READ = 'performance/read/t3-read.fasta'

T4DB = 'performance/database/t4-db.fasta'
T4READ = 'performance/read/t4-read.fasta'

T56DB = 'performance/database/polymerase-db.fasta'
T56READ = 'performance/read/polymerase-queries.fasta'

# blast and z-scores
ORDER = ['2J7W', '4K6M', '1S49', '1NB6', '3OLB', '1XR7', '1XR6', '3CDW',
         '2E9Z', '3BSO', '3UQS', '1KHV', '2CKW', '1HI0', '3AXV', '2PUS',
         '2YI9', '2R7W', '1N35', '3V81', '1MU2']

ZSCORES = {'2J7W': [60, 42.9, 22.8, 20.5, 18.1, 18.2, 18.0, 18.0, 19.2, 20.5,
                    20.9, 18.7, 17.5, 14.8, 11.1, 8.4, 9.8, 8.9, 6.5, 4.7,
                    5.4],
           '4K6M': [42.9, 60, 21.7, 17.4, 16.8, 16.6, 16.5, 16.3, 17.2, 17.5,
                    17.7, 17.9, 15.0, 10.6, 7.7, 6.6, 6.7, 9.0, 4.0, 1.6, 4.0],
           '1S49': [22.8, 21.7, 60, 27.4, 25.3, 25.1, 24.8, 25.2, 26.5, 27.1,
                    28.0, 27.4, 24.7, 4.1, 14.8, 10.7, 13.9, 10.2, 10.2, 6.3,
                    7.9],
           '1NB6': [20.5, 17.4, 27.4, 60, 21.5, 20.9, 20.7, 21.0, 21.6, 23.8,
                    25.2, 24.3, 20.6, 16.4, 14.1, 9.5, 12.9, 10.5, 7.6, 6.5,
                    7.4],
           '3OLB': [18.1, 16.8, 25.3, 21.5, 60, 52.4, 52.2, 53.1, 41.5, 32.0,
                    31.1, 32.4, 30.4, 17.2, 14.0, 12.1, 12.4, 9.7, 7.8, 5.4,
                    6.2],
           '1XR7': [18.2, 16.6, 25.1, 20.9, 52.4, 60, 56.7, 52.4, 41.3, 32.3,
                    31.5, 33.0, 30.8, 17.0, 13.5, 12.1, 12.3, 9.4, 7.3, 5.5,
                    6.6],
           '1XR6': [18.0, 16.5, 24.8, 20.7, 52.2, 56.7, 60, 53.1, 41.0, 38.1,
                    31.2, 32.9, 30.8, 16.9, 13.6, 11.9, 12.1, 8.3, 7.1, 4.9,
                    6.8],
           '3CDW': [18.0, 16.3, 25.2, 21.0, 53.1, 52.4, 53.1, 60, 41.6, 31.8,
                    31.4, 33.0, 30.9, 17.7, 14.5, 12.6, 13.0, 8.4, 7.8, 4.8,
                    6.9],
           '2E9Z': [19.2, 17.2, 26.5, 21.6, 41.5, 41.3, 41.0, 41.6, 60, 32.4,
                    32.2, 23.4, 30.8, 15.7, 13.8, 12.9, 13.5, 9.3, 8.1, 5.3,
                    6.1],
           '3BSO': [20.5, 16.5, 27.1, 23.8, 32.0, 32.3, 38.1, 31.8, 32.4, 60,
                    51.0, 39.3, 39.1, 18.5, 13.2, 13.4, 15.5, 9.4, 7.9, 5.5,
                    7.6],
           '3UQS': [20.9, 17.7, 28.0, 25.2, 31.1, 31.5, 31.2, 31.4, 32.2, 51.0,
                    60, 42.7, 39.4, 19.1, 14.4, 13.3, 14.2, 9.1, 7.9, 5.7,
                    7.9],
           '1KHV': [18.7, 17.9, 27.4, 24.3, 32.4, 33.0, 32.9, 33.0, 32.4, 39.3,
                    42.7, 60, 43.9, 17.7, 14.9, 12.6, 14.0, 10.4, 8.1, 5.7,
                    6.5],
           '1HI0': [14.8, 10.6, 4.1, 16.4, 17.2, 17.0, 16.9, 17.7, 15.7, 18.5,
                    19.1, 17.7, 14.1, 60, 12.3, 9.5, 10.7, 9.9, 8.4, 3.8, 5.5],
           '3AVX': [11.1, 7.7, 14.8, 14.1, 14.0, 13.5, 13.6, 14.5, 13.8, 13.2,
                    14.4, 14.9, 12.6, 12.3, 60, 6.0, 7.7, 7.8, 8.0, 5.8, 7.7],
           '2PUS': [8.4, 6.6, 10.7, 9.5, 12.1, 12.1, 11.9, 12.6, 12.9, 13.4,
                    13.3, 12.6, 12.9, 9.5, 6.0, 60, 42.5, 4.6, 6.5, 2.8, 3.6],
           '2YI9': [9.8, 6.7, 13.9, 12.9, 12.4, 12.3, 12.1, 13.0, 13.5, 15.5,
                    14.2, 14.0, 13.2, 10.7, 7.7, 42.5, 60, 4.6, 6.6, 2.3, 4.3],
           '2R7W': [8.9, 9.0, 10.2, 10.5, 9.7, 9.4, 8.3, 8.4, 9.3, 9.4, 9.1,
                    10.4, 8.5, 9.9, 7.8, 4.6, 4.6, 60, 15.4, 4.0, 4.6],
           '1N35': [6.5, 4.0, 10.3, 7.6, 7.8, 7.3, 7.1, 7.8, 8.1, 7.9, 7.9,
                    8.1, 8.0, 8.4, 8.0, 6.5, 6.6, 15.4, 60, 5.9, 5.1],
           '3V81': [4.7, 1.6, 6.3, 6.5, 5.4, 5.5, 4.9, 4.8, 5.3, 5.5, 5.7, 5.7,
                    4.9, 3.8, 5.8, 2.8, 2.3, 4.0, 5.9, 60, 28.5],
           '1MU2': [5.4, 4.0, 7.9, 7.4, 6.2, 6.6, 6.8, 6.9, 6.1, 7.6, 7.9, 6.5,
                    7.4, 5.5, 7.7, 3.6, 4.3, 4.6, 5.1, 28.5, 60],
           }

BITSCORES = {'2J7W': [1328.15, 899.427, 46.2098, 0, 21.9422, 0, 0, 0, 0, 0,
                      20.0162, 0, 0, 20.0162, 0, 0, 0, 0, 0, 21.557,
                      19.631],
             '4K6M': [914.835, 1901.72, 41.5874, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 20.4014, 21.557, 22.3274, 0],
             '1S49': [46.2098, 41.5874, 1268.83, 50.0618, 0, 21.1718,
                      24.2534, 0, 23.0978, 23.483, 20.4014, 25.7942, 0, 0,
                      0, 19.631, 23.8682, 24.6386, 20.7866, 0, 0],
             '1NB6': [0, 0, 50.0618, 1185.63, 0, 19.631, 0, 0, 29.6462, 0,
                      22.3274, 19.631, 0, 0, 0, 0, 0, 20.7866, 0, 0,
                      20.7866],
             '3OLB': [21.557, 0, 0, 0, 985.326, 554.288, 567.0, 719.153,
                      228.794, 55.8398, 55.8398, 69.3218, 0, 0, 0, 21.557,
                      20.4014, 21.1718, 21.557, 0, 0],
             '1XR7': [0, 0, 21.1718, 19.631, 554.288, 951.814, 840.877,
                      555.058, 218.394, 69.707, 68.1662, 52.7582, 0, 23.483,
                      0, 21.557, 19.2458, 21.557, 0, 0, 0],
             '1XR6': [0, 0, 24.2534, 0, 567.0, 840.877, 953.355, 572.392,
                      214.157, 69.3218, 63.929, 58.151, 0, 27.335, 0,
                      23.483, 21.557, 0, 0, 0, 0],
             '3CDW': [0, 0, 0, 0, 731.48, 571.622, 587.8, 978.393, 209.534,
                      51.9878, 61.2326, 58.5362, 0, 21.1718, 0, 19.2458, 0,
                      21.557, 19.2458, 0, 0],
             '2E9Z': [0, 0, 23.0978, 29.6462, 226.483, 218.394, 214.157,
                      209.149, 995.727, 61.2326, 51.9878, 47.3654, 0,
                      20.7866, 0, 0, 0, 19.631, 20.7866, 0, 23.0978],
             '3BSO': [0, 0, 23.483, 0, 55.8398, 69.707, 69.3218, 48.521,
                      61.2326, 1059.67, 619.002, 153.295, 0, 0, 0, 0, 0,
                      31.5722, 25.0238, 19.2458, 0],
             '3UQS': [20.4014, 23.8682, 20.4014, 22.3274, 55.4546, 68.5514,
                      63.929, 61.2326, 52.7582, 618.616, 1075.46, 152.525,
                      0, 0, 0, 0, 24.2534, 20.7866, 0, 22.3274, 0],
             '1KHV': [0, 19.631, 25.409, 31.187, 72.4034, 57.7658, 61.6178,
                      59.6918, 48.1358, 158.303, 162.155, 1076.62, 0, 0, 0,
                      20.0162, 0, 22.3274, 25.7942, 0, 23.8682],
             '1HI0': [20.0162, 0, 0, 0, 0, 23.483, 27.335, 0, 20.7866, 0, 0,
                      0, 0, 1385.93, 0, 0, 0, 19.631, 20.7866, 20.7866, 0],
             '3AVX': [0, 21.557, 21.557, 0, 25.7942, 26.1794, 26.1794,
                      25.7942, 24.6386, 0, 21.9422, 0, 0, 0, 0, 0, 0,
                      23.483, 0, 22.7126, 28.1054],
             '2PUS': [0, 0, 0, 0, 21.557, 21.557, 23.8682, 0, 0, 0, 0, 0,
                      0, 0, 0, 1748.02, 709.523, 20.7866, 0, 0, 0],
             '2YI9': [0, 0, 23.8682, 0, 20.4014, 0, 21.557, 0, 0, 0,
                      23.8682, 0, 0, 0, 0, 699.508, 1652.88, 25.409,
                      22.3274, 0, 20.7866],
             '2R7W': [0, 20.4014, 24.6386, 20.7866, 21.1718, 21.557, 0,
                      21.557, 0, 31.5722, 21.1718, 22.3274, 0, 0, 0,
                      20.7866, 25.409, 2261.88, 24.2534, 0, 0],
             '1N35': [0, 21.1718, 0, 0, 21.1718, 0, 0, 0, 20.7866, 25.0238,
                      0, 25.7942, 0, 20.7866, 0, 0, 22.3274, 23.8682,
                      2638.22, 21.9422, 0],
             '3V81': [21.557, 22.7126, 0, 0, 0, 0, 0, 0, 0, 0, 22.7126, 0,
                      0, 20.7866, 0, 0, 20.0162, 0, 21.9422, 1128.24,
                      690.649],
             '1MU2': [19.631, 0, 0, 20.7866, 0, 0, 0, 0, 23.0978, 0, 0,
                      23.8682, 0, 0, 0, 0, 20.7866, 0, 0, 691.419, 1137.48],
             }


def plot(x, y, read, scoreType, outputFile):
    """
    Make a scatterplot of the test results.

    @param x: a C{list} of the x coordinates (light matter score).
    @param y: a C{list} of the y coordinates (either z-score or bit score).
    @param filename: the C{str} filename where the image should be saved to.
    """
    plt.rcParams['font.family'] = 'Helvetica'
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    slope, intercept, rValue, pValue, se = stats.linregress(x, y)

    # plot
    plt.plot(x, y, 'o', markerfacecolor='blue', markeredgecolor='blue')
    if slope >= 0:
        col = 'green'
    else:
        col = 'red'
    plt.plot([0, max(x)], [intercept, slope * max(x) + intercept], '-',
             color=col)

    # labels
    ax.set_title('Read: %s, R^2: %.2f, SE: %.2f, slope: %.2f, p: %.2f' % (read,
                 rValue, se, slope, pValue))
    ax.set_ylabel(scoreType)
    ax.set_xlabel('Light matter score')

    # axes
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)

    outputDir = '/'.join(outputFile.split('/')[0:-1])
    if len(outputDir) == 0:
        outputDir = '.'
    fig.savefig('%s/%s-%s.png' % (outputDir, read, scoreType))


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

    def writeHeader(self, landmarkNames, trigNames, maxDistance, minDistance,
                    limitPerLandmark, bucketFactor):
        self.openedFile.write('Title:\nDate:\nCategory: light-matter\nTags: '
                              'light-matter, benchmarking\nSummary: '
                              'Performance and sensitivity testing\n\n'
                              '#####Database arguments:</b> '
                              'Landmarks: %s, trig points: %s, '
                              'maxDistance: %s, minDistance: %s '
                              'limitPerLandmark: %s, bucketFactor %f\n\n' %
                              (landmarkNames, trigNames,
                               maxDistance, minDistance,
                               limitPerLandmark, bucketFactor))

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
        if testName[0] == '5':
            for readName in ORDER:
                self.openedFile.write('![readName](images/%s-z.png)\n\n' % (
                                      readName))
        elif testName[0] == '6':
            for readName in ORDER:
                self.openedFile.write('![readName](images/%s-bit.png)\n\n' % (
                                      readName))

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
        '--significanceFraction', type=float, default=None,
        help='The (float) fraction of all (landmark, trig point) pairs for a '
        'scannedRead that need to fall into the same histogram bucket for '
        'that bucket to be considered a significant match with a database '
        'title.')

    parser.add_argument(
        '--verbose', default=False, action='store_true',
        help=('If True, information about each matching read will be written '
              'out.'))

    databaseSpecifier = DatabaseSpecifier(allowInMemory=False)
    databaseSpecifier.addArgsToParser(parser)

    args = parser.parse_args()

    # start writing to file
    writer = WriteMarkdownFile(args.outputFile, args.verbose)
    writer.open()
    writer.writeHeader(args.landmarkFinderNames, args.trigFinderNames,
                       args.maxDistance, args.minDistance,
                       args.limitPerLandmark, args.bucketFactor)

    # run tests
    # 1) A complete sequence must match itself:
    oneStart = time()
    print >>sys.stderr, '1) A complete sequence must find itself.'
    database = databaseSpecifier.getDatabaseFromArgs(args)
    oneResult = queryDatabase(T1DB, T1READ, database,
                              significanceFraction=args.significanceFraction)
    oneTime = time() - oneStart
    writer.writeTest('1) A complete sequence must match itself', oneResult,
                     oneTime, 1)

    # 2) Reads made from a sequence must match itself:
    twoStart = time()
    print >>sys.stderr, '2) Reads made from a sequence must match itself.'
    database = databaseSpecifier.getDatabaseFromArgs(args)
    twoResult = queryDatabase(T2DB, T2READ, database,
                              significanceFraction=args.significanceFraction)
    twoTime = time() - twoStart
    writer.writeTest('2) Reads made from a sequence must match itself',
                     twoResult, twoTime, 9)

    # 3) A sequence must find related sequences:
    threeStart = time()
    print >>sys.stderr, '3) A sequence must match related sequences.'
    database = databaseSpecifier.getDatabaseFromArgs(args)
    threeResult = queryDatabase(T3DB, T3READ, database,
                                significanceFraction=args.significanceFraction)
    threeTime = time() - threeStart
    writer.writeTest('3) A sequence must find related sequences', threeResult,
                     threeTime, 1)

    # 4) Reads made from a sequence must match related sequences:
    fourStart = time()
    print >>sys.stderr, ('4) Reads made from a sequence must match related '
                         'sequences.')
    database = databaseSpecifier.getDatabaseFromArgs(args)
    fourResult = queryDatabase(T4DB, T4READ, database,
                               significanceFraction=args.significanceFraction)
    fourTime = time() - fourStart
    fourTitle = '4) Reads made from a sequence must match related sequences'
    writer.writeTest(fourTitle, fourResult, fourTime, 7)

    # 5) The Z-scores and light matter score must correlate:
    # 6) The BLASTp bit scores and light matter scores must correlate:
    fiveSixStart = time()
    print >>sys.stderr, ('5) The Z-scores and light matter score must '
                         'correlate.')
    print >>sys.stderr, ('6) The BLASTp bit scores and light matter scores '
                         'must correlate.')
    database = databaseSpecifier.getDatabaseFromArgs(args)
    fiveSixResult = queryDatabase(
        T56DB, T56READ, database,
        significanceFraction=args.significanceFraction)
    fiveSixTime = time() - fiveSixStart
    sTitle = '6) The BLASTp bit scores and light matter scores must correlate'
    writer.writeTest('5) The Z-scores and light matter score must correlate',
                     fiveSixResult, fiveSixTime, 20)
    writer.writeTest(sTitle, fiveSixResult, fiveSixTime, 20)

    # Prepare results for plotting, and plot.
    fiveSixList = defaultdict(list)
    for read in fiveSixResult:
        for subjectName in ORDER:
            try:
                score = fiveSixResult[read][subjectName]
            except KeyError:
                score = 0
            fiveSixList[read].append(score)

    for read in fiveSixList:
        plot(fiveSixList[read], ZSCORES[read], read, 'Z-score',
             args.outputFile)
        plot(fiveSixList[read], BITSCORES[read], read, 'Bit-score',
             args.outputFile)

    # close file
    writer.close()
