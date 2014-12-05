import sys
from collections import defaultdict
import numpy as np
from scipy import stats


class ScannedReadDatabaseResult(object):
    """
    A class that holds the results from a database lookup.
    """
    def __init__(self):
        self.matches = defaultdict(lambda: defaultdict(list))
        self._finalized = False

    def __str__(self):
        if self._finalized:
            return repr(self.significant)
        else:
            raise RuntimeError('You must call finalize() before printing.')

    def addMatch(self, subjectIndex, queryId, offset):
        """
        Add a match.

        @param subjectIndex: The C{int} database index of the subject that was
            matched.
        @param queryId: The C{str} name of the query that was matched.
        @param offset: the difference in the offset of the feature in the
            subject minus the difference in the offset of the query.
        """
        self.matches[subjectIndex][queryId].append(offset)

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
        Print output.
        """
        print >>fp, self.significant
