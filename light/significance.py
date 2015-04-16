from operator import itemgetter
from warnings import warn


class Significance(object):
    """
    A class which holds functions for the calculation of significant histogram
    bins as well as the maximum significant score.

    @param histogram: A C{light.histogram} instance.
    """
    def __init__(self, histogram):
        self.histogram = histogram
        self.significantBins = []
        self.scoreGetter = itemgetter('score')

    def hashFraction(self, significanceFraction, minHashCount):
        """
        A significanceMethod which calculates the significant histogram bins
        based on which fraction of of the hashes in a query must be in a
        histogram bin.

        @param significanceFraction: The C{float} fraction of all (landmark,
            trig point) pairs for the query that need to fall into the same
            histogram bucket for that bucket to be considered a significant
            match with a database title.
        @param minHashCount: a C{int} of the total number of hashes in the
            subject or the query, depending on which one has less.
        @return: a C{list} of the significant bins.
        """
        significanceCutoff = significanceFraction * minHashCount

        # Look for bins with a significant number of elements (scaled
        # hash offset deltas).
        for binIndex, bin_ in enumerate(self.histogram.bins):
            binCount = len(bin_)
            if binCount >= significanceCutoff:
                try:
                    score = float(binCount) / minHashCount
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
                         (binCount, minHashCount), RuntimeWarning)

                self.significantBins.append({
                    'bin': bin_,
                    'index': binIndex,
                    'score': score,
                })

        return self.significantBins

    def getBestScore(self):
        """
        A method which returns the best score of a histogram. If no significant
        bins are present, None will be returned.
        """
        # Sort the significant bins by decreasing score.
        if self.significantBins:
            self.significantBins.sort(key=self.scoreGetter, reverse=True)
            return self.significantBins[0]['score']
