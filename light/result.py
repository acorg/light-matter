import sys
from collections import defaultdict
import numpy as np
from scipy import stats
from json import dumps


class ScannedReadDatabaseResult(object):
    """
    A class that holds the results from a database lookup.
    """
    def __init__(self):
        self.matches = defaultdict(lambda: defaultdict(list))
        self.matchInfo = defaultdict(lambda: defaultdict(list))
        self._finalized = False

    def __str__(self):
        if self._finalized:
            return repr(self.significant)
        else:
            raise RuntimeError('You must call finalize() before printing.')

    def addMatch(self, result):
        """
        Add a match.

        @param result: a C{dict} with information about the match.
        """
        self.matches[result['subjectIndex']][result['queryId']].append(result['combinedOffset'])
        self.matchInfo[result['subjectIndex']][result['queryId']].append({'queryOffset': result['queryOffset'],
                                                                       'subjectOffset': result['subjectOffset'],
                                                                       'subjectLength': result['length'],
                                                                       })

    def finalize(self):
        """
        Evaluates whether a subject is matched significantly by a read.
        """
        self.significant = []
        for subjectIndex in self.matches:
            for query, offsets in self.matches[subjectIndex].iteritems():
                hist, edges = np.histogram(offsets, bins=10)
                match = max(hist)
                t, p = stats.ttest_1samp(offsets, match)
                if p < 0.05:
                    self.significant.append((subjectIndex, query, match))
        self._finalized = True

    def save(self, fp=sys.stdout):
        """
        Print one line of JSON output.

        @param fp: a file pointer
        """
        alignments = []
        for subject in self.matchInfo:
            hsps = []
            query = self.matchInfo[subject]
            for hsp in self.matchInfo[subject][query]:
                hsps.append({
                    "queryOffset": hsp['queryOffset'],
                    "subjectOffset": hsp['subjectOffset'],
                })
            alignments.append({
                "hsps": hsps,
                "length": self.matchInfo[subject][query]['subjectLength'],
                "title": subject,
            })
        print >>fp, dumps({"query": query, "alignments": alignments},
                          separators=(',', ':'))
