from math import floor

# Special default value for adding to histogram to indicate that the passed
# value (to .add()) should be used as the data.
_None = object()


class Histogram(object):
    """
    Maintain a histogram.

    @param nBins: An C{int} number of bins to use in the histogram.
    """
    def __init__(self, nBins=11):
        if nBins < 1:
            raise ValueError('Number of bins must be at least one.')
        if nBins % 2 == 0:
            raise ValueError('Number of bins must be odd.')
        self.nBins = nBins
        self.max = self.min = None
        self._first = True
        self._values = []
        self._finalized = False
        self.bins = []
        for _ in range(nBins):
            self.bins.append([])

    def __getitem__(self, binIndex):
        """
        Return a particular bin, given an index.

        @param binIndex: An C{int} bin index.
        @return: A bin.
        @raise IndexError: If the bin index is out of range.
        """
        return self.bins[binIndex]

    def add(self, value, data=_None):
        """
        Add an element to the histogram.

        @param data: An object to add to the histogram.
        @param value: A numeric value associated with the object, used to
            calculate what bin the object should be placed into.
        """
        if self._finalized:
            raise RuntimeError(
                'Additional data cannot be added: histogram already finalized')
        if data is _None:
            data = value
        if self._first:
            self.max = self.min = value
            self._first = False
        else:
            if value < self.min:
                self.min = value
            elif value > self.max:
                self.max = value

        self._values.append((value, data))

    def finalize(self):
        """
        Bin all the histogram data.
        """
        if self._finalized:
            raise RuntimeError('Histogram already finalized')
        self._finalized = True
        if self._values:
            nBins = self.nBins
            min_ = self.min
            max_ = self.max
            binWidth = self.binWidth = (max_ - min_) / float(nBins)
            for value, data in self._values:
                if value == max_:
                    # The max value is put into the last bin. If we don't do
                    # this manually, we get an index error because we try to
                    # access one more bin than we have.
                    bin_ = nBins - 1
                else:
                    # Use floor to consistently round down toward
                    # numerically smaller values to get the bin number.
                    # Note that int() does not do this - it rounds towards
                    # zero!  So while floor(-3.4) = -4.0, int(-3.4) = -3.
                    bin_ = int(floor((value - min_) / binWidth))
                self.bins[bin_].append(data)
