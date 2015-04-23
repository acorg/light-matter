import numpy as np


class Always(object):
    """
    Allows every bin passed in to be significant.
    """
    def isSignificant(self, binIndex):
        """
        Determine whether a bin is significant.

        @param binIndex: The C{int} index of the bin to examine.
        @return: C{True}.
        """
        return True


class HashFraction(object):
    """
    Identify significant histogram bins based on the number of hashes in
    the bin, as compared to the theoretical maximum number, using a cutoff.

    @param histogram: A C{light.histogram} instance.
    @param minHashCount: The C{int} minimum number of hashes in the query
        and the subject.
    @param significanceFraction: The C{float} fraction of all (landmark,
        trig point) pairs for the query that need to fall into the same
        histogram bucket for that bucket to be considered a significant match
        with a database title.
    """
    def __init__(self, histogram, minHashCount, significanceFraction):
        self._histogram = histogram
        self._significanceCutoff = significanceFraction * minHashCount

    def isSignificant(self, binIndex):
        """
        Determine whether a bin is significant.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A C{bool} indicating whether the bin is significant.
        """
        binCount = len(self._histogram[binIndex])
        return binCount >= self._significanceCutoff


class MaxBinHeight(object):
    """
    Identify significant histogram bins based on the non-significant bins
    when the query is compared against itself.

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.read.AARead} instance.
    @param database: A C{light.database.Database} instance.
    """
    def __init__(self, histogram, query, database):
        self._histogram = histogram
        db = database.emptyCopy()
        subjectIndex = db.addSubject(query)
        result = db.find(query, significanceMethod='always',
                         storeFullAnalysis=True)
        bins = result.analysis[subjectIndex]['significantBins'][1:]
        binHeights = [len(h) for h in bins]
        meanBinHeight = np.mean(binHeights)
        std = np.std(binHeights)
        self._significanceCutoff = meanBinHeight + (2 * std)

    def isSignificant(self, binIndex):
        """
        Determine whether a bin is significant.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A C{bool} indicating whether the bin is significant.
        """
        binCount = len(self._histogram[binIndex])
        return binCount >= self._significanceCutoff
