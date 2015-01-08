from light.features import CombinedFeatureList


class ScannedRead(object):
    """
    Hold information about a read that has been scanned for landmarks and trig
    points.

    @param read: A C{dark.reads.Read} instance.
    """

    def __init__(self, read):
        self.read = read
        self.landmarks = []
        self.trigPoints = []

    def coveredIndices(self):
        """
        Return the indices in the read that are covered by at least one
        landmark or trig point.

        @return: The C{set} of indices that are covered by any landmark or trig
            point.
        """
        indices = set()
        for landmark in self.landmarks:
            offset = landmark.offset
            for index in xrange(offset, offset + landmark.length):
                indices.add(index)
        for trigPoint in self.trigPoints:
            offset = trigPoint.offset
            for index in xrange(offset, offset + trigPoint.length):
                indices.add(index)
        return indices

    def getPairs(self, limitPerLandmark=None, maxDistance=None,
                 minDistance=None):
        """
        Get pairs of (landmark, trig point) for use in building a search
        dictionary that can be used to identify this read.

        @param limitPerLandmark: An C{int} limit on the number of pairs to
            yield per landmark.
        @param maxDistance: The C{int} maximum distance permitted between
            yielded pairs.
        @param minDistance: The C{int} minimum distance permitted between
            yielded pairs.
        @return: A generator that yields (landmark, trig point) pairs.
        """
        if limitPerLandmark is not None and limitPerLandmark < 1:
            return

        features = CombinedFeatureList(self.landmarks, self.trigPoints)

        for landmark in self.landmarks:
            count = 0
            nearest = features.nearest(landmark.offset,
                                       maxDistance=maxDistance,
                                       minDistance=minDistance)
            while limitPerLandmark is None or count < limitPerLandmark:
                try:
                    feature = nearest.next()
                except StopIteration:
                    break
                else:
                    if feature is not landmark:
                        yield landmark, feature
                        count += 1
