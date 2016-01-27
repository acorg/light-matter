from operator import attrgetter
from functools import total_ordering

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

    def coveredOffsets(self):
        """
        Which sequence offsets occur in a feature?

        @return: A C{set} of offsets at which the feature appears.
        """
        return set(list(range(self.offset, self.offset + self.length)))


@total_ordering
class Landmark(_Feature):
    """
    Hold information about a landmark found in a sequence.

    @param name: The C{str} name of this feature.
    @param symbol: The C{str} symbol for this landmark feature.
    @param offset: The C{int} offset of the landmark in the sequence.
    @param length: The C{int} length of the landmark in the sequence.
    @param symbolDetail: Optional additional information about the symbol
        used to represent an instance of this class. For example, it
        could contain a repeat count if the landmark is made up of
        potentially many repeating sub-units. This value will be converted
        to a C{str} for use in the C{hashkey} and __lt__ functions, and so must
        have a string representation.
    """

    def __init__(self, name, symbol, offset, length, symbolDetail=''):
        _Feature.__init__(self, name, symbol, offset, length)
        self.symbolDetail = symbolDetail

    def __str__(self):
        return '%s symbol=%s offset=%d len=%d detail=%r' % (
            self.name, self.symbol, self.offset, self.length,
            self.symbolDetail)

    def __repr__(self):
        return '%s symbol=%s offset=%d len=%d detail=%r' % (
            self.name, self.symbol, self.offset, self.length,
            self.symbolDetail)

    def __eq__(self, other):
        return (self.offset == other.offset and
                self.length == other.length and
                self.symbol == other.symbol and
                self.name == other.name and
                str(self.symbolDetail) == str(other.symbolDetail))

    def __lt__(self, other):
        return (
            (self.name, self.offset, self.length, str(self.symbolDetail)) <
            (other.name, other.offset, other.length, str(other.symbolDetail)))

    def __hash__(self):
        return ('%s:%d:%d:%s' % (self.symbol, self.offset, self.length,
                                 self.symbolDetail)).__hash__()

    def hashkey(self):
        """
        Return a string for use in a hash key involving this landmark.

        @return: a C{str} of the symbol and symbol detail of the landmark.
        """
        return '%s%s' % (self.symbol, self.symbolDetail)


@total_ordering
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

    def __repr__(self):
        return '%s symbol=%s offset=%d' % (self.name, self.symbol, self.offset)

    def __eq__(self, other):
        return (self.offset == other.offset and self.symbol == other.symbol and
                self.name == other.name)

    def __lt__(self, other):
        return (self.name, self.offset) < (other.name, other.offset)

    def __hash__(self):
        return ('%s:%d' % (self.symbol, self.offset)).__hash__()

    def hashkey(self):
        """
        Return a string for use in a hash key involving this trig point.

        @return: a C{str} of the symbol for the trig point.
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
        self.landmarks = landmarks
        self.trigPoints = trigPoints

    def nearest(self, offset, maxDistance=None, minDistance=None):
        """
        Find the features nearest to the given offset. The first features
        returned will be landmarks, and then trig points.

        @param offset: The C{int} offset whose nearest features are wanted.
        @param maxDistance: The C{int} maximum distance a feature may be from
            the given offset to be considered.
        @param minDistance: The C{int} minimum distance between a feature and
            the given offset.
        @return: A generator that yields nearby features.
        """
        key = attrgetter('offset')

        # Sort the landmarks separately first and return them. We will sort
        # the trig points only if we need to. Remember that we are returning
        # a generator and our caller may not call next() enough times to
        # trigger the second (trig point) sort.

        features = SortedCollection(self.landmarks, key=key)
        for feature in self._nearest(offset, features, maxDistance,
                                     minDistance):
            yield feature

        features = SortedCollection(self.trigPoints, key=key)
        for feature in self._nearest(offset, features, maxDistance,
                                     minDistance):
            yield feature

    def _nearest(self, offset, features, maxDistance, minDistance):
        """
        Helper function for finding the features in a given sorted collection
        nearest to a given offset.

        @param offset: The C{int} offset whose nearest features are wanted.
        @param features: An instance of C{SortedCollection}.
        @param maxDistance: The C{int} maximum distance a feature may be from
            the given offset to be considered.
        @param minDistance: The C{int} minimum distance between a feature and
            the given offset.
        @return: A generator that yields nearby features.
        """

        nFeatures = len(features)

        # Set right to be the index of the first feature with offset
        # greater than or equal to the wanted offset.
        try:
            rightFeature = features.find_ge(offset)
        except ValueError:
            right = nFeatures
        else:
            right = features.index(rightFeature)

        left = right - 1

        while True:
            if left >= 0:
                leftDelta = abs(features[left].offset - offset)
            else:
                leftDelta = None
            if right < nFeatures:
                rightDelta = abs(features[right].offset - offset)
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
                            yield features[right]
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
                            yield features[left]
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
                                yield features[left]
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
                                yield features[right]
                                right += 1
                            else:
                                # the offset is smaller than the minDistance,
                                # don't yield, but continue.
                                right += 1
                        else:
                            # The smallest available delta is too large.
                            return
