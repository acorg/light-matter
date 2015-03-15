class Histogram(object):
    """
    Maintain a histogram.

    @param nBins: An C{int} number of bins to use in the histogram.
    """
    def __init__(self, nBins=10):
        if nBins < 1:
            raise ValueError('Number of bins must be at least one.')
        self.max = self.min = None
        self._first = True
        self._values = []
        self._finalized = False
        self.bins = []
        for _ in xrange(nBins):
            self.bins.append([])

    def add(self, data, value):
        """
        Add an element to the histogram.

        @param data: An object to add to the histogram.
        @param value: A numeric value associated with the object, used to
            calculate what bin the object should be place into.
        """
        assert not self._finalized
        if self._first:
            self.max = self.min = value
            self._first = False
        else:
            if value < self.min:
                self.min = value
            elif value > self.max:
                self.max = value

        self._values.append((data, value))

    def finalizeHistogram(self):
        """
        Bin all the histogram data.
        """
        self._finalized = True
        # Bin the data.
        nBins = len(self.bins)
        binWidth = self.binWidth = (self.max - self.min) / float(nBins)
        min_ = self.min
        max_ = self.max
        for data, value in self._values:
            if value == max_:
                # The max value is put into the last bin. If we don't do
                # this manually, we get an index error because we try to
                # access one more bin than we have.
                bin_ = nBins - 1
            else:
                bin_ = int((value - min_) / binWidth)
            self.bins[bin_].append(data)
