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

    def addMatch(self, info, subjectIndex):
        """
        Add a match.

        @param info: a C{dict} with information about the landmark and
            trigpoint names, offsets, landmark length and distance between
            landmark and trigpoint.
        @param subjectIndex: a C{int} index of the subject in the database.
        """
        if subjectIndex in self.matches:
            self.matches[subjectIndex]['info'].append(info)
        else:
            self.matches[subjectIndex] = {
                'info': [info],
            }

    def finalize(self, aboveMeanThreshold):
        """
        Evaluates whether a subject is matched significantly by a read.

        @param aboveMeanThreshold: A numeric amount by which the maximum delta
            count in a bucket must exceed the mean bucket count for that
            maximum bucket count to be considered significant.
        """
        for subjectIndex in self.matches:
            offsets = [info['offsets']['subjectOffset'] -
                       info['offsets']['readOffset']
                       for info in self.matches[subjectIndex]['info']]
            hist, edges = np.histogram(offsets)
            mean = np.mean(hist)
            match = max(hist)
            if match >= mean + aboveMeanThreshold:
                self.matches[subjectIndex]['matchScore'] = match
                self.significant[subjectIndex] = self.matches[subjectIndex]
        self._finalized = True

    def save(self, fp=sys.stdout):
        """
        Print one line of JSON output.

        @param fp: a file pointer.
        """
        print 'noooooot'
        alignments = []
        for subjectIndex in self.significant:
            hsps = self.significant[subjectIndex]['info']
            matchScore = self.significant[subjectIndex]['matchScore']
            alignments.append({
                'hsps': hsps,
                'matchScore': matchScore,
                'subjectIndex': subjectIndex,
            })
        print >>fp, dumps({
                          'alignments': alignments,
                          'query': self.read.id,
                          'querySequence': self.read.sequence,
                          },
                          separators=(',', ':'))
