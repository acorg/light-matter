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
