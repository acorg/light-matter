import sys
from json import dumps
from collections import defaultdict
from warnings import warn

from light.distance import scale
from light.histogram import Histogram


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
    @param hashCount: The C{int} number of hashes that the query has (including
        those that were not found in the database).
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
    """
    def __init__(self, scannedQuery, matches, hashCount, significanceFraction,
                 database, nonMatchingHashes=None,
                 storeFullAnalysis=False):
        self.scannedQuery = scannedQuery
        self.matches = matches  # Only saved on self for testing.
        self.hashCount = hashCount
        self.significanceFraction = significanceFraction
        self.database = database
        self.nonMatchingHashes = nonMatchingHashes
        self._storeFullAnalysis = storeFullAnalysis
        self.analysis = defaultdict(dict)
        significanceCutoff = significanceFraction * hashCount
        distanceBase = database.distanceBase
        queryLen = len(scannedQuery.read.sequence)

        # Go through all the subejcts that were matched at all, and put the
        # match offset deltas into bins so we can decide which of the
        # matches is significant.
        for subjectIndex in matches:
            subjectLen = len(database.subjectInfo[subjectIndex][1])
            # Use a histogram to bin scaled (landmark, trigPoint) offset
            # deltas.
            maxLen = max([queryLen, subjectLen])
            nBins = scale(maxLen, distanceBase)
            histogram = Histogram(nBins)
            add = histogram.add

            for match in matches[subjectIndex]:
                landmark = match['landmark']
                trigPoint = match['trigPoint']

                for queryLandmarkOffset, queryTrigPointOffset in zip(
                        match['queryLandmarkOffsets'],
                        match['queryTrigPointOffsets']):
                    for subjectOffset in match['subjectLandmarkOffsets']:
                        delta = subjectOffset - queryLandmarkOffset
                        add({
                            'landmark': landmark,
                            'queryLandmarkOffset': queryLandmarkOffset,
                            'queryTrigPointOffset': queryTrigPointOffset,
                            'subjectLandmarkOffset': subjectOffset,
                            'trigPoint': trigPoint,
                            },
                            scale(delta, distanceBase))

            histogram.finalizeHistogram()
            significantBins = []
            bestScore = None

            # Look for bins with a significant number of elements (deltas).
            for binIndex, bin_ in enumerate(histogram.bins):
                count = len(bin_)
                if count > significanceCutoff:
                    # We don't have to worry about hashCount being zero here.
                    # Because 'matches' has some hashes in it, hashCount
                    # cannot be zero by this point.
                    score = count / float(hashCount)

                    # It is possible that a bin has more deltas in it than
                    # the number of hashes in the query. This can occur
                    # when the distanceBase used to scale distances results
                    # in two or more different distances being put into the
                    # same bin. In this case, the distance binning is
                    # arguably being too aggessive. For now, issue a
                    # warning and set the score to 1.0
                    if score > 1.0:
                        score = 1.0
                        warn('Bin contains %d deltas for a query that has '
                             'only %d hashes. The database distanceBase is '
                             'causing different distances to be scaled into '
                             'the same bin, which might not be a good thing.' %
                             (count, hashCount), RuntimeWarning)

                    significantBins.append({
                        'bin': bin_,
                        'score': score,
                    })
                    if bestScore is None or score > bestScore:
                        bestScore = score

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

    def significant(self):
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
        for subjectIndex in self.significant():
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
               printFeatures=True, printHistograms=False,
               queryDescription='Query title'):
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

        result = [
            'Overall matches: %d' % len(subjectIndices),
            'Significant matches: %d' % len(list(self.significant())),
            'Query hash count: %d' % self.hashCount,
            'Significance fraction: %f' % self.significanceFraction,
            'Significance cutoff: %f' % (self.significanceFraction *
                                         self.hashCount),
        ]

        if subjectIndices:
            result.append('Matched subjects:')

        for subjectCount, subjectIndex in enumerate(subjectIndices, start=1):
            analysis = self.analysis[subjectIndex]
            title, sequence = self.database.subjectInfo[subjectIndex]

            # Get a list of the significant bins, sorted by decreasing score.
            significantBins = sorted(
                analysis['significantBins'], reverse=True,
                key=lambda bin_: bin_['score'])

            result.extend([
                '  Subject %d title: %s' % (subjectCount, title),
                '    Index in database: %d' % subjectIndex,
                '    Best HSP score: %s' % analysis['bestScore'],
            ])

            if printSequences:
                result.append('    Sequence: %s' % sequence)

            result.append('    HSP count: %d' % len(significantBins))

            for hspCount, bin_ in enumerate(significantBins, start=1):
                pairCount = len(bin_['bin'])
                result.append('      HSP %d has %d matching (landmark, '
                              'trigpoint) pair%s and score: %f' % (
                                  hspCount, pairCount,
                                  '' if pairCount == 1 else 's',
                                  bin_['score']))
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
                binCounts = [len(b) for b in histogram.bins]
                result.extend([
                    '    Histogram:',
                    '      Number of bins: %d' % len(histogram.bins),
                    '      Max bin count: %r' % max(binCounts),
                    '      Counts: %r' % binCounts,
                    '      Max (scaled) offset delta: %d' % histogram.max,
                    '      Min (scaled) offset delta: %d' % histogram.min,
                ])

        print >>fp, '\n'.join(result)
