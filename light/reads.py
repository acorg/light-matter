from light.features import CombinedFeatureList, TrigPoint


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

        When considering each landmark, other landmarks to the left (with
        smaller offset) of the current landmark are ignored.  This is to
        reduce the redundancy that results if both (landmark1, landmark2)
        and (landmark2, landmark1) pairs are yielded.

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
            landmarkOffset = landmark.offset
            nearest = features.nearest(landmarkOffset,
                                       maxDistance=maxDistance,
                                       minDistance=minDistance)

            while limitPerLandmark is None or count < limitPerLandmark:
                try:
                    feature = nearest.next()
                except StopIteration:
                    break
                else:
                    # Yield any landmark or trig point on the right (with
                    # greater offset) of the current landmark, and all trig
                    # points (whether left or right).
                    if (feature.offset > landmarkOffset or
                            isinstance(feature, TrigPoint)):
                        yield landmark, feature
                        count += 1
