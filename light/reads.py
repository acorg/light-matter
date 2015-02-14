from light.features import CombinedFeatureList, Landmark


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
        smaller offset) of the current landmark will only be yielded after
        other possible pairs have been yielded. This is to reduce the
        redundancy that results if both (landmark1, landmark2) and
        (landmark2, landmark1) pairs are yielded.

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
            landmarksOnLeft = []
            nearest = features.nearest(landmarkOffset,
                                       maxDistance=maxDistance,
                                       minDistance=minDistance)

            # Take features from the nearest iterator, but do not yield any
            # features that are landmarks that are to the left (i.e., with
            # smaller offset) of the current landmark. Instead, save these
            # and only yield them after we have yielded everything else
            # (landmarks to the right, trig points on either side).
            while limitPerLandmark is None or count < limitPerLandmark:
                try:
                    feature = nearest.next()
                except StopIteration:
                    break
                else:
                    if feature is not landmark:
                        if (feature.offset < landmarkOffset and
                                isinstance(feature, Landmark)):
                            landmarksOnLeft.append(feature)
                        else:
                            yield landmark, feature
                            count += 1

            # Now yield landmarks to the left of the current landmark, if any.
            for landmarkOnLeft in landmarksOnLeft:
                if limitPerLandmark is None or count < limitPerLandmark:
                    yield landmark, landmarkOnLeft
                    count += 1
                else:
                    break
