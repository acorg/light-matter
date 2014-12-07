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
            off = offsets['subjectOffset'] - offsets['readOffset']
            self.matches[subjectIndex]['offsets'].append(off)
            self.matches[subjectIndex]['absoluteOffsets'].append(offsets)
        else:
            off = offsets['subjectOffset'] - offsets['readOffset']
            self.matches[subjectIndex] = {
                'subjectLength': subjectLength,
                'offsets': [off],
                'absoluteOffsets': [offsets],
            }

    def finalize(self):
        """
        Evaluates whether a subject is matched significantly by a read.
        """
        self.significant = []
        for subjectIndex in self.matches:
            offsets = self.matches[subjectIndex]['offsets']
            hist, edges = np.histogram(offsets, bins=10)
            match = max(hist)
            t, p = stats.ttest_1samp(offsets, match)
            if p < 0.05:
                self.significant.append((subjectIndex, self.read.id, match))
        self._finalized = True

    def save(self, fp=sys.stdout):
        """
        Print one line of JSON output.

        @param fp: a file pointer.
        """
        alignments = []
        for subjectIndex in self.matches:
            hsps = []
            for absoluteOffsets in self.matches[subjectIndex]['absoluteOffsets']:
                hsps.append({
                    "readOffset": absoluteOffsets['readOffset'],
                    "subjectOffset": absoluteOffsets['subjectOffset'],
                })
            alignments.append({
                "hsps": hsps,
                "subjectLength": self.matches[subjectIndex]['subjectLength'],
                "subjectIndex": subjectIndex,
            })
        print >>fp, dumps({"query": self.read.id, "alignments": alignments},
                          separators=(',', ':'))
