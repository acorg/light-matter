import sys
import numpy as np
from json import dumps


class Result(object):
    """
    A class that holds the results from a database lookup.

    @param read: a C{dark.read.AARead} instance.
    @param database: A C{light.database.Database} instance.
    """
    def __init__(self, read, database):
        self.matches = {}
        self.read = read
        self._database = database
        self._finalized = False
        self.significant = {}

    def __str__(self):
        if self._finalized:
            return repr(self.significant)
        else:
            raise RuntimeError('You must call finalize() before printing.')

    def addMatch(self, offsets, subjectIndex):
        """
        Add a match.

        @param offsets: a C{dict} with information about the match.
        @param subjectIndex: a C{int} index of the subject in the database.
        """
        if subjectIndex in self.matches:
            self.matches[subjectIndex]['offsets'].append(offsets)
        else:
            self.matches[subjectIndex] = {
                'offsets': [offsets],
            }

    def finalize(self):
        """
        Evaluates whether a subject is matched significantly by a read.
        """
        for subjectIndex in self.matches:
            offsets = [offsets['subjectOffset'] - offsets['readOffset']
                       for offsets in self.matches[subjectIndex]['offsets']]
            hist, edges = np.histogram(offsets)
            mean = np.mean(hist)
            match = max(hist)
            if mean + 15 < match:
                self.matches[subjectIndex]['matchScore'] = match
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
            subjectTitle, subjectSequence = self._database.subjectInfo[
                subjectIndex]
            matchScore = self.significant[subjectIndex]['matchScore']
            alignments.append({
                'hsps': hsps,
                'matchScore': matchScore,
                'subjectLength': len(subjectSequence),
                'subjectIndex': subjectIndex,
                'subjectTitle': subjectTitle,
            })
        print >>fp, dumps({
                          'alignments': alignments,
                          'query': self.read.id,
                          'querySequence': self.read.sequence,
                          },
                          separators=(',', ':'))
