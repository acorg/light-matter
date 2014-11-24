class Landmark(object):
    """
    Hold information about a landmark found in a sequence.

    @param symbol: The C{str} symbol for this landmark feature.
    @param offset: The C{int} offset of the landmark in the sequence.
    @param repeatCount: The C{int} number of times the landmark pattern was
        found in the sequence at this offset.
    """

    def __init__(self, symbol, offset, repeatCount=1):
        self.symbol = symbol
        self.offset = offset
        self.repeatCount = repeatCount

    def __str__(self):
        return '%s%d at %d' % (self.symbol, self.repeatCount, self.offset)

    def __eq__(self, other):
        return (self.symbol == other.symbol and self.offset == other.offset and
                self.repeatCount == other.repeatCount)
