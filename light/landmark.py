class Landmark(object):
    """
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
