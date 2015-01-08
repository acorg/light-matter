from operator import attrgetter

from light.sortedCollection import SortedCollection


class _Feature(object):
    """
    Hold information about a landmark or trig point found in a sequence.

    @param name: The C{str} name of this feature.
    @param symbol: The C{str} symbol for this feature.
    @param offset: The C{int} offset of the feature in the sequence.
    """

    def __init__(self, name, symbol, offset, length):
        self.name = name
        self.symbol = symbol
        self.offset = offset
        self.length = length


class Landmark(_Feature):
    """
    Hold information about a landmark found in a sequence.

    @param name: The C{str} name of this feature.
    @param symbol: The C{str} symbol for this landmark feature.
    @param offset: The C{int} offset of the landmark in the sequence.
    @param length: The C{int} length of the landmark in the sequence.
    @param repeatCount: The C{int} number of times the landmark pattern was
        found in the sequence at this offset.
    """

    def __init__(self, name, symbol, offset, length, repeatCount=1):
        _Feature.__init__(self, name, symbol, offset, length)
        self.repeatCount = repeatCount

    def __str__(self):
        return '%s symbol=%s%d offset=%d len=%d' % (
            self.name, self.symbol, self.repeatCount, self.offset, self.length)

    def __eq__(self, other):
        return (self.symbol == other.symbol and
                self.offset == other.offset and
                self.length == other.length and
                self.repeatCount == other.repeatCount)

    def hashkey(self):
        """
        Return a string suitable for use as a hash key for this landmark.

        @return: a C{str} of the symbol of the respective landmark.
        """
        return '%s%d' % (self.symbol, self.repeatCount)


class TrigPoint(_Feature):
    """
    Hold information about a trigonometric point found in a sequence.

    @param name: The C{str} name of this feature.
    @param symbol: The C{str} symbol for this trig point.
    @param offset: The C{int} offset of the trig point in the sequence.
    """

    def __init__(self, name, symbol, offset):
        _Feature.__init__(self, name, symbol, offset, 1)

    def __str__(self):
        return '%s symbol=%s offset=%d' % (self.name, self.symbol, self.offset)

    def __eq__(self, other):
        return (self.name == other.name and
                self.symbol == other.symbol and
                self.offset == other.offset)

    def hashkey(self):
        """
        Return a string suitable for use as a hash key for this trig point.

        @return: a C{str} of the symbol of the respective trig point.
        """
        return self.symbol


class CombinedFeatureList(object):
    """
    Hold a sorted (by offset) collection of _Feature instances (landmarks
    and trig points), and provide a method for finding nearby features.

    @param landmarks: An iterable of L{Landmark} instances.
    @param trigPoints: An iterable of L{TrigPoint} instances.
    """
    def __init__(self, landmarks, trigPoints):
        self._features = SortedCollection(landmarks + trigPoints,
                                          key=attrgetter('offset'))

    def nearest(self, offset, maxDistance=None, minDistance=None):
        """
        Find the features nearest to the given offset.

        @param offset: The C{int} offset whose nearest features are wanted.
        @param maxDistance: The C{int} maximum distance a feature may be from
            the given offset to be considered.
        @param minDistance: The C{int} minimum distance between a feature and
            the given offset.
        @return: A generator that yields nearby features.
        """

        nFeatures = len(self._features)

        # Set right to be the index of the first feature with offset
        # greater than or equal to the wanted offset.
        try:
            rightFeature = self._features.find_ge(offset)
        except ValueError:
            right = nFeatures
        else:
            right = self._features.index(rightFeature)

        left = right - 1

        while True:
            if left >= 0:
                leftDelta = abs(self._features[left].offset - offset)
            else:
                leftDelta = None
            if right < nFeatures:
                rightDelta = abs(self._features[right].offset - offset)
            else:
                rightDelta = None

            if leftDelta is None:
                if rightDelta is None:
                    # We ran out of neighboring features on left and right.
                    return
                else:
                    # We only have a right neighboring feature.
                    if maxDistance is None or rightDelta <= maxDistance:
                        if minDistance is None or rightDelta >= minDistance:
                            yield self._features[right]
                            right += 1
                        else:
                            # the offset is smaller than the minDistance,
                            # don't yield, but continue.
                            right += 1
                    else:
                        # The next delta is too large. We're done.
                        return

            else:
                if rightDelta is None:
                    # We only have a left neighboring feature.
                    if maxDistance is None or leftDelta <= maxDistance:
                        if minDistance is None or leftDelta >= minDistance:
                            yield self._features[left]
                            left -= 1
                        else:
                            # the offset is smaller than the minDistance,
                            # don't yield, but continue.
                            left -= 1
                    else:
                        # The next delta is too large. We're done.
                        return

                else:
                    # Deltas on both left and right are available.
                    if leftDelta < rightDelta:
                        if maxDistance is None or leftDelta <= maxDistance:
                            if minDistance is None or leftDelta >= minDistance:
                                yield self._features[left]
                                left -= 1
                            else:
                                # the offset is smaller than the minDistance,
                                # don't yield, but continue.
                                left -= 1
                        else:
                            # The smallest available delta is too large.
                            return

                    else:
                        if maxDistance is None or rightDelta <= maxDistance:
                            if (minDistance is None or
                                    rightDelta >= minDistance):
                                yield self._features[right]
                                right += 1
                            else:
                                # the offset is smaller than the minDistance,
                                # don't yield, but continue.
                                right += 1
                        else:
                            # The smallest available delta is too large.
                            return
