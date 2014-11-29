class _Feature(object):
    """
    Hold information about a landmark or trig point found in a sequence.

    @param name: The C{str} name of this feature.
    @param symbol: The C{str} symbol for this feature.
    @param offset: The C{int} offset of the feature in the sequence.
    @param repeatCount: The C{int} number of times the feature's pattern was
        found in the sequence at this offset.
    """

    def __init__(self, name, symbol, offset, length, repeatCount=1):
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
