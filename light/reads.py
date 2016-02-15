from light.features import CombinedFeatureList, TrigPoint
from light.string import MultilineString


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
            for index in range(offset, offset + landmark.length):
                indices.add(index)
        for trigPoint in self.trigPoints:
            offset = trigPoint.offset
            # Trig points are always length one, but use trigPoint.length
            # in case that ever changes.
            for index in range(offset, offset + trigPoint.length):
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
                    feature = next(nearest)
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

    def print_(self, printSequence=False, printFeatures=False,
               description='Read', margin='', result=None):
        """
        Print the details of a scanned read.

        @param fp: A file pointer to print to.
        @param verbose: If C{True}, print details of landmark and trig
            point matches.
        @param printSequence: If C{True}, print the sequence.
        @param printFeatures: If C{True}, print details of landmark and trig
            point features.
        @param description: A C{str} description to print before the scanned
            read id. This allows us to be specific about what a read is, e.g.,
            a query or a subject.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @param result: A C{MultilineString} instance, or C{None} if a new
            C{MultilineString} should be created.
        @return: If C{result} was C{None}, return a C{str} representation of
            the scanned read, else C{None}.
        """
        read = self.read

        if result is None:
            result = MultilineString(margin=margin)
            returnNone = False
        else:
            returnNone = True

        append = result.append
        coveredIndices = len(self.coveredIndices())

        append('%s: %s' % (description, read.id))
        result.indent()
        if printSequence:
            append('Sequence: %s' % read.sequence)
        append('Length: %d' % len(read.sequence))
        append('Covered indices: %d (%.2f%%)' % (
            coveredIndices,
            coveredIndices / float(len(read.sequence)) * 100.0))

        # Print read landmarks and trig points.
        append('Landmark count %d, trig point count %d' % (
            len(self.landmarks), len(self.trigPoints)))
        if printFeatures:
            result.indent()
            for landmark in self.landmarks:
                append(str(landmark))
            for trigPoint in self.trigPoints:
                append(str(trigPoint))
            result.outdent()

        result.outdent()

        if not returnNone:
            return str(result)
