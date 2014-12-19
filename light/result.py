import sys
import numpy as np
from json import dumps


class Result(object):
    """
    Hold the result of a database lookup on a read.

    @param read: A C{dark.read.AARead} instance.
    @param matches: A C{dict} of matches. Keys are C{int} subject indices.
        Each value is a C{list} of C{dicts}, with each C{dict} containing the
        following keys: 'distance', 'landmarkLength', 'landmarkName',
        'readOffset', 'subjectOffset', and 'trigPointName'.
    @param aboveMeanThreshold: A numeric amount by which the maximum count
        across all buckets must exceed the mean bucket count for the
        maximum bucket count to be considered significant.
    """
    def __init__(self, read, matches, aboveMeanThreshold):
        self.matches = matches
        self.read = read
        self.significant = set()
        self.scores = {}
        for subjectIndex in matches:
            offsets = [match['subjectOffset'] - match['readOffset']
                       for match in matches[subjectIndex]]
            hist = np.histogram(offsets)[0]
            mean = np.mean(hist)
            match = max(hist)
            self.scores[subjectIndex] = match
            if match >= mean + aboveMeanThreshold:
                self.significant.add(subjectIndex)

    def __str__(self):
        return repr(self.significant)

    def save(self, fp=sys.stdout):
        """
        Print a line of JSON with the significant results for this read.

        @param fp: a file pointer.
        """
        alignments = []
        for subjectIndex in self.significant:
            alignments.append({
                'hsps': [{
                    'matchInfo': self.matches[subjectIndex],
                    'matchScore': self.scores[subjectIndex],
                }],
                'subjectIndex': subjectIndex,
            })
        print >>fp, dumps(
            {
                'alignments': alignments,
                'queryId': self.read.id,
                'querySequence': self.read.sequence,
            },
            separators=(',', ':'))
