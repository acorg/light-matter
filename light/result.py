import sys
import numpy as np
from scipy import stats
from json import dumps


class ScannedReadDatabaseResult(object):
    """
    A class that holds the results from a database lookup.

    @param read: a C{dark.read.AARead} instance.
    """
    def __init__(self, read):
        self.matches = {}
        self.read = read
        self._finalized = False
        self.significant = {}

    def __str__(self):
        if self._finalized:
            return repr(self.significant)
        else:
            raise RuntimeError('You must call finalize() before printing.')

    def addMatch(self, offsets, subjectIndex, subjectLength):
        """
        Add a match.

        @param offsets: a C{dict} with information about the match.
        @param subjectIndex: a C{int} index of the subject in the database.
        @param subjectLength: a C{int} length of the subject.
        """
        if subjectIndex in self.matches:
            self.matches[subjectIndex]['offsets'].append(offsets)
        else:
            self.matches[subjectIndex] = {
                'subjectLength': subjectLength,
                'offsets': [offsets],
            }

    def finalize(self):
        """
        Evaluates whether a subject is matched significantly by a read.
        """
        for subjectIndex in self.matches:
            offsets = [offsets['subjectOffset'] - offsets['readOffset']
                       for offsets in self.matches[subjectIndex]['offsets']]
            hist, edges = np.histogram(offsets, bins=10)
            match = max(hist)
            t, p = stats.ttest_1samp(offsets, match)
            if p < 0.05:
                self.significant[subjectIndex] = self.matches[subjectIndex]
        self._finalized = True

    def save(self, fp=sys.stdout):
        """
        Print one line of JSON output.

        @param fp: a file pointer.
        """
        alignments = []
        for subjectIndex in self.significant:
            hsps = self.significant[subjectIndex]['offsets']
            alignments.append({
                'hsps': hsps,
                'subjectLength': self.significant[subjectIndex]['subjectLength'],
                'subjectIndex': subjectIndex,
            })
        print >>fp, dumps({'query': self.read.id, 'alignments': alignments},
                          separators=(',', ':'))
