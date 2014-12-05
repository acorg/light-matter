import sys
from collections import defaultdict
import numpy as np
from scipy import stats


class ScannedReadDatabaseResult(object):
    """
    A class that holds the results from a database lookup.

    @param subjectId: a C{str} name of the subject that was matched.
    @param queryId: a C{str} name of the query that was matched.
    @param offset: the difference in the offset of the feature in the subject
        minus the difference in the offset of the query
    """
    def __init__(self):
        self.matches = defaultdict(lambda: defaultdict(list))
        self.finalized = False

    def __str__(self):
        if self.finalized:
            return repr(self.significant)
        else:
            pass
            # 'not finalized'

    def addMatch(self, subjectId, queryId, offset):
        """
        Add a match.

        @param subjectId: a C{str} name of the subject that was matched.
        @param queryId: a C{str} name of the query that was matched.
        @param offset: the difference in the offset of the feature in the
            subject minus the difference in the offset of the query.
        """
        self.matches[subjectId][queryId].append(offset)

    def finalize(self):
        """
        Evaluates whether a subject is matched significantly by a read.
        """
        self.significant = []
        for subject in self.matches:
            for query, offsets in self.matches[subject].iteritems():
                hist, edges = np.histogram(offsets, bins=10)
                match = max(hist)
                t, p = stats.ttest_1samp(offsets, match)
                if p < 0.05:
                    self.significant.append((subject, query, match))
        self.finalized = True

    def save(self, fp=sys.stdout):
        """
        Print output.
        """
        print >>fp, self
