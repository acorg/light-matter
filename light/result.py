import sys
from copy import deepcopy
from json import dumps
from collections import defaultdict
from math import log10, ceil
from operator import itemgetter
from warnings import warn

from light.distance import scale
from light.histogram import Histogram
from light.significance import HashFraction, MaxBinHight


class Result(object):
    """
    Hold the result of a database find() for a read.

    @param scannedQuery: A C{dark.read.ScannedRead} instance.
    @param matches: A C{dict} of match information. Keys are the C{int}
        database subject indices of subjects that the scanned read had some
        hashes in common with. Each value is a C{list} of C{dicts}. Each
        of those C{dict}s corresponds to one hash that occurs in both the
        query and the subject, and contains keys: 'landmark', 'readOffsets',
        'subjectOffsets', and 'trigPoint'.
    @param queryHashCount: The C{int} number of hashes that the query has
        (including those that were not found in the database).
    @param significanceMethod: The name of the method used to calculate
        which histogram bins are considered significant.
    @param significanceFraction: The C{float} fraction of all (landmark,
        trig point) pairs for the query that need to fall into the same
        histogram bucket for that bucket to be considered a significant match
        with a database title.
    @param database: A C{light.database.Database} instance.
    @param nonMatchingHashes: A C{dict} whose keys are hashes that were found
        in C{scannedQuery} but that did not match any subject. Each value is a
        C{dict} with keys: 'landmark', 'readOffsets', and 'trigPoint'.
    @param storeFullAnalysis: A C{bool}. If C{True} the full significance
        analysis of each matched subject will be stored.
    @raises ValueError: If a non-existing significanceMethod is specified.
    """
    def __init__(self, scannedQuery, matches, queryHashCount,
                 significanceMethod, significanceFraction, database,
                 nonMatchingHashes=None, storeFullAnalysis=False):
        self.scannedQuery = scannedQuery
        self.matches = matches  # Only saved on self for testing.
        self.queryHashCount = queryHashCount
        self.significanceMethod = significanceMethod
        self.significanceFraction = significanceFraction
        self.database = database
        self.nonMatchingHashes = nonMatchingHashes
        self._storeFullAnalysis = storeFullAnalysis
        self.analysis = defaultdict(dict)
        distanceBase = database.distanceBase
        queryLen = len(scannedQuery.read.sequence)
        scoreGetter = itemgetter('score')

        # Go through all the subjects that were matched at all, and put the
        # match offset deltas into bins so we can decide which (if any) of
        # the matches is significant.
        for subjectIndex in matches:
            subject = database.getSubject(subjectIndex)
            # Use a histogram to bin scaled (landmark, trigPoint) offset
            # deltas.
            subjectLen = len(subject)
            maxLen = max([queryLen, subjectLen])
            nBins = scale(maxLen, distanceBase)
            # Make sure the number of bins is odd, else Histogram() will raise.
            nBins |= 0x1
            histogram = Histogram(nBins)
            add = histogram.add
            # To ensure the set of query/subject offset deltas is the same
            # no matter which of the sequences is the query and which is
            # the subject, we negate all deltas if the subject sequence
            # sorts first.  This is just a way of canonicalizing the set of
            # deltas.  If we don't canonicalize, we get sets of deltas with
            # opposite signs, like {-4, -2, 6} and {-6, 2, 4} depending on
            # which sequence is the subject and which the query. This
            # occasionally leads to hard-to-debug and awkward-to-fix
            # differences in the histogram binning at bin boundaries due to
            # tiny floating point differences. The simple solution is to
            # canonicalize the deltas based on an arbitrary consistent
            # difference between the subject and query.
            negateDeltas = subject.sequence < scannedQuery.read.sequence

            for match in matches[subjectIndex]:
                landmark = match['landmark']
                trigPoint = match['trigPoint']

                for queryLandmarkOffset, queryTrigPointOffset in zip(
                        match['queryLandmarkOffsets'],
                        match['queryTrigPointOffsets']):
                    for subjectOffset in match['subjectLandmarkOffsets']:
                        delta = subjectOffset - queryLandmarkOffset
                        if negateDeltas:
                            delta = -delta
                        add(scale(delta, distanceBase),
                            {
                                'landmark': landmark,
                                'queryLandmarkOffset': queryLandmarkOffset,
                                'queryTrigPointOffset': queryTrigPointOffset,
                                'subjectLandmarkOffset': subjectOffset,
                                'trigPoint': trigPoint,
                            })

            histogram.finalize()

            minHashCount = min(queryHashCount, subject.hashCount)

            if significanceMethod == 'hashFraction':
                significance = HashFraction(
                    histogram, minHashCount, significanceFraction)
            elif significanceMethod == 'maxBinHight':
                significance = MaxBinHight(histogram, scannedQuery.read,
                                           database)
            else:
                raise ValueError('Unknown significance method %r' %
                                 significanceMethod)
            # Look for bins with a significant number of elements (each
            # element is a scaled hash offset delta).
            significantBins = []
            for binIndex, bin_ in enumerate(histogram.bins):
                if significance.isSignificant(binIndex):
                    binCount = len(bin_)
                    try:
                        score = float(binCount) / minHashCount
                    except ZeroDivisionError:
                        score = 0.0

                    # It is possible that a bin has more deltas in it than
                    # the number of hashes in the query / subject. This can
                    # occur when the distanceBase used to scale distances
                    # results in multiple different distances being put
                    # into the same bin. In this case, the distance binning
                    # is perhaps too aggessive. For now, issue a warning
                    # and set the score to 1.0
                    if score > 1.0:
                        score = 1.0
                        warn('Bin contains %d deltas for a query/subject pair '
                             'with a minimum hash count of only %d. The '
                             'database distanceBase is causing different '
                             'distances to be scaled into the same bin, which '
                             'might not be a good thing.' %
                             (binCount, minHashCount), RuntimeWarning)

                    significantBins.append({
                        'bin': bin_,
                        'index': binIndex,
                        'score': score,
                    })

            if significantBins:
                significantBins.sort(key=scoreGetter, reverse=True)
                bestScore = significantBins[0]['score']
            else:
                bestScore = None

            if storeFullAnalysis:
                self.analysis[subjectIndex] = {
                    'histogram': histogram,
                    'bestScore': bestScore,
                    'significantBins': significantBins,
                }
            elif significantBins:
                self.analysis[subjectIndex] = {
                    'bestScore': bestScore,
                    'significantBins': significantBins,
                }

    def significantSubjects(self):
        """
        Which subject matches were significant?

        @return: A generator that yields the subject indices of the
            significant matched subjects.
        """
        if self._storeFullAnalysis:
            return (subjectIndex for subjectIndex, analysis in
                    self.analysis.iteritems()
                    if analysis['significantBins'])
        else:
            # When the full analysis isn't being stored, all keys (i.e.,
            # subject indices) in self.analysis were significant.
            return self.analysis.iterkeys()

    def save(self, fp=sys.stdout):
        """
        Print a line of JSON with the significant subject matches for this
        read.

        @param fp: a file pointer.
        @return: The C{fp} we were passed (this is useful in testing).
        """
        alignments = []
        for subjectIndex in self.significantSubjects():
            analysis = self.analysis[subjectIndex]
            hsps = []

            for significantBin in analysis['significantBins']:
                # Each significant bin for a subject is what BLAST calls an
                # HSP: a significant way in which the query can be aligned
                # to the subject.
                hspInfo = []
                for binItem in significantBin['bin']:
                    hspInfo.append({
                        'landmark': binItem['landmark'].name,
                        'landmarkLength': binItem['landmark'].length,
                        'queryLandmarkOffset': binItem['queryLandmarkOffset'],
                        'queryTrigPointOffset':
                        binItem['queryTrigPointOffset'],
                        'subjectLandmarkOffset':
                        binItem['subjectLandmarkOffset'],
                        'trigPoint': binItem['trigPoint'].name,
                    })

                hsps.append({
                    'hspInfo': hspInfo,
                    'score': significantBin['score'],
                })

            alignments.append({
                'hsps': hsps,
                'matchScore': analysis['bestScore'],
                'subjectIndex': subjectIndex,
            })

        read = self.scannedQuery.read
        print >>fp, dumps(
            {
                'alignments': alignments,
                'queryId': read.id,
                'querySequence': read.sequence,
            },
            separators=(',', ':'))

        return fp

    def print_(self, fp=sys.stdout, printQuery=True, printSequences=False,
               printFeatures=False, printHistograms=False,
               queryDescription='Query title', sortHSPsByScore=True):
        """
        Print a result in a human-readable format. If self._storeFullAnalysis
        is True, full information about all matched subjects (i.e., including
        matches that were not significant) will be printed. If not, only basic
        information about significant matches will appear.

        @param fp: A file pointer to print to.
        @param printQuery: If C{True}, also print details of the query.
        @param printSequences: If C{True}, also print query and subject
            sequences.
        @param printFeatures: If C{True}, print details of landmark and trig
            point features.
        @param printHistograms: If C{True}, print details of histograms.
        @param queryDescription: A C{str} description to print before the query
            (when printQuery is C{True}.
        @param sortHSPsByScore: If C{True}, HSPs for a subject should be
            printed in order of decreasing score. If C{False}, print sorted by
            histogram bin number.
        """
        if printQuery:
            self.scannedQuery.print_(fp=fp, printSequence=printSequences,
                                     printFeatures=printFeatures,
                                     description=queryDescription)

        # Sort matched subjects (if any) in order of decreasing score so we
        # can print them in a useful order.
        #
        # The following sorted() call will fail (with TypeError) under
        # Python 3 because bestScore is None when there are no significant
        # matches (which can happen when self._storeFullAnalysis is True).
        subjectIndices = sorted(
            self.analysis.iterkeys(), reverse=True,
            key=lambda index: self.analysis[index]['bestScore'])

        if not sortHSPsByScore:
            indexGetter = itemgetter('index')

        result = [
            'Overall matches: %d' % len(subjectIndices),
            'Significant matches: %d' % len(list(self.significantSubjects())),
            'Query hash count: %d' % self.queryHashCount,
            'Significance fraction: %f' % self.significanceFraction,
        ]

        if subjectIndices:
            result.append('Matched subjects:')

        for subjectCount, subjectIndex in enumerate(subjectIndices, start=1):
            analysis = self.analysis[subjectIndex]
            subject = self.database.getSubject(subjectIndex)
            minHashCount = min(self.queryHashCount, subject.hashCount)
            significantBins = analysis['significantBins']

            result.extend([
                '  Subject %d:' % subjectCount,
                '    Title: %s' % subject.id,
                '    Best HSP score: %s' % analysis['bestScore'],
            ])

            if printSequences:
                result.append('    Sequence: %s' % subject.sequence)

            result.extend([
                '    Index in database: %d' % subjectIndex,
                '    Subject hash count: %s' % subject.hashCount,
                '    Subject/query min hash count: %s' % minHashCount,
                '    Significance cutoff: %f' % (self.significanceFraction *
                                                 minHashCount),
                '    Number of HSPs: %d' % len(significantBins),
            ])

            if not sortHSPsByScore:
                significantBins = deepcopy(significantBins)
                significantBins.sort(key=indexGetter)

            for hspCount, bin_ in enumerate(significantBins, start=1):
                binCount = len(bin_['bin'])
                result.append(
                    '      HSP %d (bin %d): %d matching hash%s, score %f' %
                    (hspCount, bin_['index'], binCount,
                     '' if binCount == 1 else 'es', bin_['score']))

                if printFeatures:
                    for binItem in bin_['bin']:
                        result.extend([
                            '        Landmark %s subjectOffset=%d' % (
                                binItem['landmark'],
                                binItem['subjectLandmarkOffset']),
                            '        Trig point %s' % binItem['trigPoint'],
                        ])

            if printHistograms and self._storeFullAnalysis:
                histogram = analysis['histogram']
                significantBinIndices = set([bin_['index']
                                             for bin_ in significantBins])
                maxCount = max(len(bin_) for bin_ in histogram.bins)

                result.extend([
                    '    Histogram:',
                    '      Number of bins: %d' % len(histogram.bins),
                    '      Bin width: %.10f' % histogram.binWidth,
                    '      Max bin count: %r' % maxCount,
                    '      Max (scaled) offset delta: %d' % histogram.max,
                    '      Min (scaled) offset delta: %d' % histogram.min,
                ])

                result.extend([
                ])

                # Calculate column widths for displaying ranges neatly.
                maxAbsoluteValue = max(
                    [-histogram.min, histogram.max, -histogram.max])
                if maxAbsoluteValue == 0:
                    # All printed range values will be '+0.0', of length 4.
                    rangeWidth = 4
                else:
                    # Add 3 because we have the sign, a decimal point, and
                    # one digit of precision.
                    rangeWidth = 3 + int(ceil(log10(maxAbsoluteValue)))
                rangeSeparator = ' to '
                rangeColumnWidth = 2 * rangeWidth + len(rangeSeparator)

                first = True
                for binIndex, bin_ in enumerate(histogram.bins):
                    binCount = len(bin_)
                    if binCount:
                        if first:
                            result.extend([
                                '      Non-empty bins:',
                                '        %s %s %*s %s' % (
                                    'Index', 'Count', rangeColumnWidth,
                                    'Range', 'Significant'),
                            ])
                            first = False
                        binLow = histogram.min + binIndex * histogram.binWidth
                        # The 5, 5, 11 below are the lengths of 'Index',
                        # 'Range', and 'Significant'.
                        result.append('        %5d %5d %+*.1f%s%+*.1f %11s' % (
                            binIndex, binCount,
                            rangeWidth, binLow,
                            rangeSeparator,
                            rangeWidth, binLow + histogram.binWidth,
                            'Yes' if binIndex in significantBinIndices
                            else ''))
                if first:
                    result.append('All bins were empty.')

        print >>fp, '\n'.join(result)
