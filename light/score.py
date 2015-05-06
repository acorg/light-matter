from warnings import warn


class MinHashesScore(object):
    """
    Calculates the score for a bin based on the minimum hashes of either
    subject or query.

    @param histogram: A C{light.histogram} instance.
    @param minHashCount: The C{int} minimum number of hashes in the query
        and the subject.
    """
    def __init__(self, histogram, minHashCount):
        self._histogram = histogram
        self.minHashCount = minHashCount

    def calculateScore(self, binIndex):
        """
        Calculates the score.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A C{float} of the score of that bin.
        """
        binCount = len(self._histogram[binIndex])
        try:
            score = float(binCount) / self.minHashCount
        except ZeroDivisionError:
            score = 0.0

        # It is possible that a bin has more deltas in it than
        # the number of hashes in the query / subject. This can
        # occur when the distanceBase used to scale distances
        # results in multiple different distances being put
        # into the same bin. In this case, the distance binning
        # is perhaps too aggessive. For now, issue a warning
        # and set the score to 1.0
        if score > 1.0:
            score = 1.0
            warn('Bin contains %d deltas for a query/subject pair '
                 'with a minimum hash count of only %d. The '
                 'database distanceBase is causing different '
                 'distances to be scaled into the same bin, which '
                 'might not be a good thing.' %
                 (binCount, self.minHashCount), RuntimeWarning)

        return score


ALL_SCORE_CLASSES = (MinHashesScore,)
