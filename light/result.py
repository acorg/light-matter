import sys
from json import dumps
from collections import defaultdict

from light.histogram import Histogram


class Result(object):
    """
    Hold the result of a database look-up on a read.

    @param scannedRead: A C{dark.read.ScannedRead} instance.
    @param matches: A C{dict} of matches. Keys are C{int} subject indices.
        Each value is a C{list} of C{dicts}, with each C{dict} containing the
        following keys: 'landmark', 'subjectOffsets', and 'trigPoint'.
    @param hashCount: The C{int} number of hashes that a scannedRead has
        (including hashes that were not found in the database).
    @param significanceFraction: The C{float} fraction of all (landmark,
        trig point) pairs for a scannedRead that need to fall into the
        same histogram bucket for that bucket to be considered a
        significant match with a database title.
    @param bucketFactor: A C{int} factor by which the distance between
        landmark and trig point is divided, to influence sensitivity.
    @param nonMatchingHashes: A C{set} of hashes in scannedRead that does not
        match a subject.
    @param storeFullAnalysis: A C{bool}. If C{True} the full significance
        analysis of each matched subject will be stored.
    """
    def __init__(self, scannedRead, matches, hashCount, significanceFraction,
                 bucketFactor, nonMatchingHashes=None,
                 storeFullAnalysis=False):
        self.scannedRead = scannedRead
        self.matches = matches
        self._storeFullAnalysis = storeFullAnalysis
        self.analysis = defaultdict(dict)
        significanceCutoff = significanceFraction * hashCount
        for subjectIndex in matches:
            maxLen = max([len(scannedRead.read.sequence),
                          matches[subjectIndex][0]['subjectLength']])
            nBins = maxLen // bucketFactor
            histogram = Histogram(nBins)

            for match in matches[subjectIndex]:
                readOffset = match['landmark'].offset
                for subjectOffset in match['subjectOffsets']:
                    delta = subjectOffset - readOffset
                    histogram.add(match, delta)

            histogram.finalizeHistogram()
            significant = [bin_ for bin_ in histogram.bins
                           if len(bin_) > significanceCutoff]

            if significant:
                pairs = set()
                for bin_ in significant:
                    for match in bin_:
                        pairs.add((match['landmark'], match['trigPoint']))
                score = len(pairs) / float(hashCount)
            else:
                score = None

            if storeFullAnalysis:
                self.analysis[subjectIndex] = {
                    'hashCount': hashCount,
                    'histogram': histogram,
                    'score': score,
                    'significanceCutoff': significanceCutoff,
                    'significantBinCount': len(significant),
                }
            elif significant:
                self.analysis[subjectIndex] = {
                    'score': score,
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
                    if analysis['significantBinCount'])
        else:
            # When the full analysis isn't being stored, all keys in
            # self.analysis are significant.
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
            matches = []
            for match in self.matches[subjectIndex]:
                matches.append({
                    'landmarkLength': match['landmark'].length,
                    'landmarkName': match['landmark'].name,
                    'readOffset': match['landmark'].offset,
                    'subjectLength': match['subjectLength'],
                    'subjectOffsets': match['subjectOffsets'],
                    'trigPointName': match['trigPoint'].name,
                })
            alignments.append({
                'matchInfo': matches,
                'matchScore': self.analysis[subjectIndex]['score'],
                'subjectIndex': subjectIndex,
            })

        read = self.scannedRead.read
        print >>fp, dumps(
            {
                'alignments': alignments,
                'queryId': read.id,
                'querySequence': read.sequence,
            },
            separators=(',', ':'))

        return fp
