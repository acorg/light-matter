import sys
from json import dumps
from collections import defaultdict

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
        self.matches = matches
        self.hashCount = hashCount
        self.significanceFraction = significanceFraction
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
                    score = count / float(hashCount)
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

    def print_(self, database, fp=sys.stdout, printQuery=True, verbose=False):
        """
        Print a result in a human-readable format. If self._storeFullAnalysis
        is True, full information about all matched subjects (i.e., including
        matches that were not significant) will be printed. If not, only basic
        information about significant matches will appear.

        @param database: A C{light.database.Database} instance.
        @param fp: A file pointer to print to.
        @param printQuery: If C{True}, also print details of the query.
        @param verbose: If C{True}, print details of landmark and trig
            point matches.
        """
        if printQuery:
            self.scannedQuery.print_(fp, verbose)

        print >>fp, 'Significant matches: %d' % len(list(self.significant()))
        print >>fp, 'Overall matches: %d' % len(self.matches)
        print >>fp, 'Hash count: %d' % self.hashCount
        print >>fp, 'Significance fraction: %f' % self.significanceFraction
        print >>fp, 'Significance cutoff: %f' % (self.significanceFraction *
                                                 self.hashCount)
        if self.matches:
            print >>fp, 'Subject matches:'
            # Print matched subjects in order of decreasing score.
            subjectIndices = sorted(
                self.analysis.iterkeys(), reverse=True,
                key=lambda index: self.analysis[index]['bestScore'])

            for subjectIndex in subjectIndices:
                title, sequence = database.subjectInfo[subjectIndex]
                analysis = self.analysis[subjectIndex]
                print >>fp, '  Title: %s' % title
                print >>fp, '    Score: %s' % analysis['bestScore']
                print >>fp, '    Sequence: %s' % sequence
                print >>fp, '    Database subject index: %d' % subjectIndex

                if self._storeFullAnalysis:
                    histogram = analysis['histogram']
                    binCounts = [len(b) for b in histogram.bins]
                    print >>fp, '    Histogram'
                    print >>fp, '      Number of bins: %d' % (
                        len(histogram.bins))
                    print >>fp, '      Significant bin count: %s' % (
                        len(analysis['significantBins']))
                    print >>fp, '      Max bin count: %r' % max(binCounts)
                    print >>fp, '      Counts: %r' % binCounts
                    print >>fp, '      Max (scaled) offset delta: %d' % (
                        histogram.max)
                    print >>fp, '      Min (scaled) offset delta: %d' % (
                        histogram.min)
