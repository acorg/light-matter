import sys
import numpy as np
from json import dumps
from collections import defaultdict


class Result(object):
    """
    Hold the result of a database look-up on a read.

    @param read: A C{dark.read.AARead} instance.
    @param matches: A C{dict} of matches. Keys are C{int} subject indices.
        Each value is a C{list} of C{dicts}, with each C{dict} containing the
        following keys: 'landmarkLength', 'landmarkName',
        'readOffset', 'subjectOffsets', and 'trigPointName'.
    @param hashCount: The C{int} number of hashes that a scannedRead has
        (including hashes that were not found in the database).
    @param significanceFraction: The C{float} fraction of all (landmark,
        trig point) pairs for a scannedRead that need to fall into the
        same histogram bucket for that bucket to be considered a
        significant match with a database title.
    @param bucketFactor: A C{int} factor by which the distance between
        landmark and trig point is divided, to influence sensitivity.
    @param storeAnalysis: A C{bool}. If C{True} the intermediate significance
        analysis of each matched subject will be stored. Else it is discarded.
    """
    def __init__(self, read, matches, hashCount, significanceFraction,
                 bucketFactor, storeAnalysis=False):
        self.read = read
        self.matches = matches
        self.analysis = defaultdict(dict)
        for subjectIndex in matches:
            offsets = [subjectOffset - match['readOffset']
                       for match in matches[subjectIndex]
                       for subjectOffset in match['subjectOffsets']]
            maxLen = max([len(read.sequence),
                          matches[subjectIndex][0]['subjectLength']])
            bins = maxLen // bucketFactor
            histogram, histogramBuckets = np.histogram(offsets, bins=bins)
            maxCount = np.max(histogram)
            maxCountFraction = maxCount / float(hashCount)
            significant = (maxCountFraction >= significanceFraction)
            self.analysis[subjectIndex] = {
                'score': maxCount,
                'significant': significant,
            }

            if storeAnalysis:
                self.analysis[subjectIndex].update({
                    'hashCount': hashCount,
                    'histogram': histogram,
                    'histogramBuckets': histogramBuckets,
                    'maxCount': maxCount,
                    'offsets': offsets,
                })

    def significant(self):
        """
        Which subject matches were significant?

        @return: A generator that yields the subject indices of the
            significant matched subjects.
        """
        return (subjectIndex for subjectIndex, analysis in
                self.analysis.iteritems() if analysis['significant'])

    def save(self, fp=sys.stdout):
        """
        Print a line of JSON with the significant results for this read.

        @param fp: a file pointer.
        @return: The C{fp} we were passed (this is useful in testing).
        """
        alignments = []
        for subjectIndex in self.significant():
            alignments.append({
                'matchInfo': self.matches[subjectIndex],
                'matchScore': self.analysis[subjectIndex]['score'],
                'subjectIndex': subjectIndex,
            })

        print >>fp, dumps(
            {
                'alignments': alignments,
                'queryId': self.read.id,
                'querySequence': self.read.sequence,
            },
            separators=(',', ':'))

        return fp
