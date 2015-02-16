from unittest import TestCase

from light.histogram import Histogram


class TestHistogram(TestCase):
    """
    Tests for the light.histogram.Histogram class
    """

    def testNoData(self):
        """
        If no data is given to the Histogram class, all bins must be empty.
        """
        h = Histogram()
        for bin_ in h.bins:
            self.assertEqual(0, len(bin_))

    def testNoDataMaxAndMin(self):
        """
        If no data is given to the Histogram class, the max and min attributes
        of the histogram must be None.
        """
        h = Histogram()
        self.assertIs(None, h.max)
        self.assertIs(None, h.min)

    def testDefaultNumberOfBins(self):
        """
        If no bin count is given to the Histogram class, it must make the
        expected number of bins.
        """
        h = Histogram()
        self.assertEqual(10, len(h.bins))

    def testNonDefaultNumberOfBins(self):
        """
        If a bin count is given to the Histogram class, it must make the
        expected number of bins.
        """
        h = Histogram(5)
        self.assertEqual(5, len(h.bins))

    def testOneElementOneBin(self):
        """
        If a histogram is created with just one element and one bin, the
        element must be placed in the bin.
        """
        h = Histogram(1)
        h.add(3, 3)
        h.finalizeHistogram()
        self.assertEqual([[3]], h.bins)

    def testOneElementMaxMin(self):
        """
        If a histogram is created with just one element, the max and min
        should be set to that value.
        """
        h = Histogram()
        h.add(3, 3)
        h.finalizeHistogram()
        self.assertEqual(3, h.max)
        self.assertEqual(3, h.min)

    def testTwoElementsMaxMin(self):
        """
        If a histogram is created with two elements, the max and min
        should be set to the correct values.
        """
        h = Histogram()
        h.add(3, 3)
        h.add(3, 4)
        h.finalizeHistogram()
        self.assertEqual(4, h.max)
        self.assertEqual(3, h.min)

    def testThreeElementsMaxMin(self):
        """
        If a histogram is created with three elements, the max and min
        should be set to the correct values.
        """
        h = Histogram()
        h.add(3, 3)
        h.add(3, 4)
        h.add(3, 5)
        h.finalizeHistogram()
        self.assertEqual(5, h.max)
        self.assertEqual(3, h.min)

    def testElementIsStoredInBin(self):
        """
        If a histogram is created with just one element and one bin, the
        exact element that was passed must be placed in the bin.
        """
        element = object()
        h = Histogram(1)
        h.add(element, 3)
        h.finalizeHistogram()
        self.assertIs(element, h.bins[0][0])

    def testTenElementsInTwoBins(self):
        """
        If a histogram is created with ten elements placed into 2 bins, the
        bins must contain the expected values.
        """
        h = Histogram(2)
        for i in range(10):
            h.add(i, i)
        h.finalizeHistogram()
        self.assertEqual([[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]], h.bins)
