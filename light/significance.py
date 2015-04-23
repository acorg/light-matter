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
