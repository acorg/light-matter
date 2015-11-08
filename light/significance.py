import numpy as np


class Always(object):
    """
    Every bin is significant.
    """
    def __init__(self):
        self.analysis = {
            'significanceMethod': self.__class__.__name__,
        }

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
        self.significanceCutoff = significanceFraction * minHashCount
        self.analysis = {
            'significanceMethod': self.__class__.__name__,
            'significanceCutoff': self.significanceCutoff,
        }

    def isSignificant(self, binIndex):
        """
        Determine whether a bin is significant.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A C{bool} indicating whether the bin is significant.
        """
        binCount = len(self._histogram[binIndex])
        return binCount >= self.significanceCutoff


class MaxBinHeight(object):
    """
    Identify significant histogram bins based on the maximum non-significant
    bins when the query is compared against itself.

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.read.AARead} instance.
    @param database: A C{light.database.Database} instance.
    """
    def __init__(self, histogram, query, database):
        self._histogram = histogram
        # A top-level import of Database would be circular.
        from light.database import Database
        db = Database(database.params)
        _, subjectIndex, _ = db.addSubject(query)
        from light.parameters import FindParameters
        findParams = FindParameters(significanceMethod='Always')
        result = db.find(query, findParams, storeFullAnalysis=True)
        bins = result.analysis[subjectIndex]['histogram'].bins
        # The highest-scoring bin is ignored.
        binHeights = sorted([len(h) for h in bins], reverse=True)[1:]
        self.significanceCutoff = binHeights[0]
        self.analysis = {
            'significanceMethod': self.__class__.__name__,
            'significanceCutoff': self.significanceCutoff,
        }

    def isSignificant(self, binIndex):
        """
        Determine whether a bin is significant.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A C{bool} indicating whether the bin is significant.
        """
        binCount = len(self._histogram[binIndex])
        return binCount >= self.significanceCutoff


class MeanBinHeight(object):
    """
    Identify significant histogram bins based on the mean non-significant bins
    when the query is compared against itself.

    @param histogram: A C{light.histogram} instance.
    @param query: A C{dark.read.AARead} instance.
    @param database: A C{light.database.Database} instance.
    """
    def __init__(self, histogram, query, database):
        self._histogram = histogram
        # A top-level import of Database would be circular.
        from light.database import Database
        db = Database(database.params)
        _, subjectIndex, _ = db.addSubject(query)
        from light.parameters import FindParameters
        findParams = FindParameters(significanceMethod='Always')
        result = db.find(query, findParams, storeFullAnalysis=True)
        bins = result.analysis[subjectIndex]['histogram'].bins
        # The highest-scoring bin is ignored.
        binHeights = sorted([len(h) for h in bins], reverse=True)[1:]
        self.meanBinHeight = np.mean(binHeights)
        self.std = np.std(binHeights)
        self.significanceCutoff = self.meanBinHeight + (2 * self.std)
        self.analysis = {
            'significanceMethod': self.__class__.__name__,
            'significanceCutoff': self.significanceCutoff,
            'standardDeviation': self.std,
            'meanBinHeight': self.meanBinHeight,
        }

    def isSignificant(self, binIndex):
        """
        Determine whether a bin is significant.

        @param binIndex: The C{int} index of the bin to examine.
        @return: A C{bool} indicating whether the bin is significant.
        """
        binCount = len(self._histogram[binIndex])
        return binCount >= self.significanceCutoff


ALL_SIGNIFICANCE_CLASSES = (Always, HashFraction, MaxBinHeight, MeanBinHeight)
