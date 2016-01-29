import six
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

    def testFinalizeWithNoData(self):
        """
        If finalize is called and the histogram has no data, all
        bins must have zero counts.
        """
        h = Histogram(5)
        h.finalize()
        self.assertEqual([0, 0, 0, 0, 0], [len(bin_) for bin_ in h.bins])

    def testZeroBins(self):
        """
        If a bin count of zero is given to the Histogram class, a ValueError
        must be raised.
        """
        error = '^Number of bins must be at least one\.$'
        six.assertRaisesRegex(self, ValueError, error, Histogram, 0)

    def testEvenNumberOfBins(self):
        """
        If an even bin count is given to the Histogram class, a ValueError
        must be raised.
        """
        error = '^Number of bins must be odd\.$'
        six.assertRaisesRegex(self, ValueError, error, Histogram, 6)

    def testNegativeBins(self):
        """
        If a negative bin count is given to the Histogram class, a ValueError
        must be raised.
        """
        error = '^Number of bins must be at least one\.$'
        six.assertRaisesRegex(self, ValueError, error, Histogram, -1)

    def testAddDataAfterFinalized(self):
        """
        If an attempt is made to add to a histogram that has been finalized,
        a RuntimeError must be raised.
        """
        h = Histogram()
        error = ('^Additional data cannot be added: histogram already '
                 'finalized$')
        h.add(3)
        h.finalize()
        six.assertRaisesRegex(self, RuntimeError, error, h.add, 3)

    def testRepeatedFinalize(self):
        """
        If finalize is called a second time, a RuntimeError must be raised.
        """
        h = Histogram()
        error = ('^Histogram already finalized$')
        h.add(3)
        h.finalize()
        six.assertRaisesRegex(self, RuntimeError, error, h.finalize)

    def testDefaultNumberOfBins(self):
        """
        If no bin count is given to the Histogram class, it must make the
        expected default number of bins.
        """
        h = Histogram()
        self.assertEqual(11, h.nBins)

    def testNonDefaultNumberOfBins(self):
        """
        The given number of bins must be stored on the histogram.
        """
        h = Histogram(17)
        self.assertEqual(17, h.nBins)

    def testActualNumberOfBins(self):
        """
        If a bin count is given to the Histogram class, it must make the
        expected number of bins.
        """
        h = Histogram(5)
        self.assertEqual(5, len(h.bins))

    def testElementIsStoredInBin(self):
        """
        If a histogram is created with just one element and one bin, the
        exact element that was passed must be placed in the bin.
        """
        element = object()
        h = Histogram(1)
        h.add(3, element)
        h.finalize()
        self.assertIs(element, h.bins[0][0])

    def testNoDataValue(self):
        """
        If an element with no associated datum is added to a histogram,
        the value that is passed must be stored in the bin.
        """
        h = Histogram(1)
        h.add(3)
        h.finalize()
        self.assertEqual([[3]], h.bins)

    def testOneElementBinWidth(self):
        """
        If a histogram is created with just one element, the bin width must be
        zero.
        """
        h = Histogram()
        h.add(3)
        h.finalize()
        self.assertEqual(0.0, h.binWidth)

    def testOneElementMaxMin(self):
        """
        If a histogram is created with just one element, the max and min
        should be set to that value.
        """
        h = Histogram()
        h.add(3)
        h.finalize()
        self.assertEqual(3, h.max)
        self.assertEqual(3, h.min)

    def testTwoElementsBinWidth(self):
        """
        If a histogram with 5 buckets is created with two elements that differ
        by 1.0, the bin width should be set to the correct value of 0.2.
        """
        h = Histogram(5)
        h.add(3)
        h.add(4)
        h.finalize()
        self.assertEqual(0.2, h.binWidth)

    def testTwoElementsMaxMin(self):
        """
        If a histogram is created with two elements, the max and min
        should be set to the correct values.
        """
        h = Histogram()
        h.add(3)
        h.add(4)
        h.finalize()
        self.assertEqual(4, h.max)
        self.assertEqual(3, h.min)

    def testThreeElementsMaxMin(self):
        """
        If a histogram is created with three elements, the max and min
        should be set to the correct values.
        """
        h = Histogram()
        h.add(3)
        h.add(4)
        h.add(5)
        h.finalize()
        self.assertEqual(5, h.max)
        self.assertEqual(3, h.min)

    def testNineElementsInThreeBins(self):
        """
        If a histogram is created with 9 elements placed into 2 bins, the
        bins must contain the expected values.
        """
        h = Histogram(3)
        list(map(h.add, range(9)))
        h.finalize()
        self.assertEqual([[0, 1, 2], [3, 4, 5], [6, 7, 8]], h.bins)

    def testTenElementsInThreeBinsBinWidth(self):
        """
        If a histogram is created with 10 elements (0-9) placed into 3 bins,
        the bin width must be 3.0.
        """
        h = Histogram(3)
        list(map(h.add, range(10)))
        h.finalize()
        self.assertEqual(3.0, h.binWidth)

    def testFiveBinsMinusTwoPointFiveToPlusTwoPointFiveIntermediates(self):
        """
        If a histogram is created with 5 bins and a data range of -2.5 to +2.5
        items that are added between histogram boundaries must be placed in
        the expected bins.
        """
        for (value, expectedCounts) in ((-2, [1, 0, 0, 0, 0]),
                                        (-1, [0, 1, 0, 0, 0]),
                                        (+0, [0, 0, 1, 0, 0]),
                                        (+1, [0, 0, 0, 1, 0]),
                                        (+2, [0, 0, 0, 0, 1])):
            h = Histogram(5)
            h.add(-2.5)  # Set min value.
            h.add(2.5)  # Set max value.
            h.add(value)
            h.finalize()
            counts = [len(bin_) for bin_ in h.bins]
            # Subract 1 from the first and last bin counts, to adjust for the
            # -2.5 and 2.5 boundary values we added manually.
            counts[0] -= 1
            counts[-1] -= 1
            self.assertEqual(expectedCounts, counts)

    def testFiveBinsMinusTwoPointFiveToPlusTwoPointFiveBoundaries(self):
        """
        If a histogram is created with 5 bins and a data range of -2.5 to +2.5
        items that are added must be placed in the expected bins when the
        values fall on the histogram boundaries.
        """
        for (value, expectedCounts) in ((-2.5, [1, 0, 0, 0, 0]),
                                        (-1.5, [0, 1, 0, 0, 0]),
                                        (-0.5, [0, 0, 1, 0, 0]),
                                        (+0.5, [0, 0, 0, 1, 0]),
                                        (+1.5, [0, 0, 0, 0, 1]),
                                        (+2.5, [0, 0, 0, 0, 1])):
            h = Histogram(5)
            h.add(-2.5)  # Set min value.
            h.add(2.5)  # Set max value.
            h.add(value)
            h.finalize()
            counts = [len(bin_) for bin_ in h.bins]
            # Subract 1 from the first and last bin counts, to adjust for the
            # -2.5 and 2.5 boundary values we added manually.
            counts[0] -= 1
            counts[-1] -= 1
            self.assertEqual(expectedCounts, counts)

    def _checkPositiveNegative(self, nBins, values):
        """
        When a set of values is put into a histogram, the bin counts that
        result must be the same (just with the order reversed) as those that
        result from a histogram made with the same set of values but with
        opposite sign.

        @param nBins: The C{int} number of bins to use in the histogram.
        @param values: A C{list} of values to insert into the histogram.
        """
        # Make a histogram of the values and get all the bin counts.
        h1 = Histogram(nBins)
        for value in values:
            h1.add(value)
        h1.finalize()
        counts1 = [len(bin_) for bin_ in h1.bins]

        # Make a histogram of the negative values and get all the bin counts.
        h2 = Histogram(nBins)
        for value in [-x for x in values]:
            h2.add(value)
        h2.finalize()
        counts2 = [len(bin_) for bin_ in h2.bins]
        counts2.reverse()

        # Prepare a useful error message, in case there are any differences.
        differences = ['Counts differ']
        for i in range(len(counts1)):
            if counts1[i] != counts2[i]:
                h1Low = h1.min + i * h1.binWidth
                h1High = h1Low + h1.binWidth
                h2Low = h2.min + i * h2.binWidth
                h2High = h2Low + h2.binWidth
                differences.append(
                    '  bin %d (h1 bin range: %.7f to %.7f, h2 bin range: '
                    '%.7f to %.7f): count %d != count %d' % (
                        i, h1Low, h1High, h2Low, h2High,
                        counts1[i], counts2[i]))

        # Bin counts must be the same.
        self.assertEqual(counts1, counts2, '\n'.join(differences))

    def testBinCountsPostiveNegative(self):
        """
        See docstring for self._checkPositiveNegative (above).

        This example is taken from work on
        https://github.com/acorg/light-matter/issues/235
        """
        values = [
            -262, -262, -262, -262, -262, -247, -247, -247, -247, -247,
            -247, -247, -247, -233, -233, -233, -233, -228, -228, -228, -228,
            -223, -223, -223, -223, -223, -223, -223, -223, -223, -216, -216,
            -216, -216, -216, -216, -216, -216, -211, -211, -211, -211, -211,
            -211, -211, -211, -211, -207, -207, -207, -207, -207, -207, -207,
            -207, -207, -207, -207, -203, -203, -203, -192, -192, -192, -192,
            -192, -192, -192, -192, -192, -178, -178, -178, -178, -178, -178,
            -178, -178, -178, -178, -164, -164, -164, -164, -164, -164, -164,
            -164, -164, -164, -164, -164, -164, -162, -162, -162, -162, -162,
            -158, -158, -158, -158, -158, -158, -158, -158, -158, -152, -152,
            -152, -152, -152, -152, -152, -152, -152, -152, -148, -148, -148,
            -148, -148, -148, -134, -134, -134, -134, -134, -134, -134, -134,
            -134, -134, -134, -134, -133, -133, -133, -133, -133, -129, -129,
            -129, -129, -129, -128, -128, -128, -128, -128, -128, -126, -126,
            -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -121,
            -121, -121, -121, -121, -121, -121, -121, -121, -121, -121, -109,
            -109, -109, -109, -109, -109, -109, -109, -109, -109, -109, -109,
            -109, -109, -109, -109, -103, -103, -103, -103, -103, -102, -102,
            -102, -102, -102, -102, -90, -90, -90, -90, -90, -90, -90, -90,
            -87, -87, -87, -87, -87, -87, -87, -87, -87, -87, -87, -87, -80,
            -80, -80, -80, -77, -77, -77, -76, -76, -76, -76, -65, -65, -65,
            -65, -65, -65, -65, -64, -64, -64, -64, -64, -64, -62, -62, -62,
            -62, -62, -62, -62, -62, -62, -62, -62, -62, -59, -59, -59, -59,
            -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59,
            -59, -59, -59, -59, -59, -59, -59, -59, -56, -56, -56, -56, -56,
            -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56,
            -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56,
            -56, -56, -56, -56, -56, -56, -56, -55, -55, -55, -55, -55, -55,
            -55, -55, -53, -53, -53, -53, -53, -53, -52, -52, -52, -52, -52,
            -47, -47, -47, -47, -47, -47, -47, -47, -47, -44, -44, -44, -44,
            -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44,
            -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44,
            -44, -44, -44, -44, -44, -44, -44, -40, -40, -40, -40, -40, -40,
            -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40,
            -40, -39, -39, -39, -39, -39, -39, -39, -31, -31, -30, -30, -30,
            -30, -30, -30, -30, -30, -30, -29, -29, -29, -29, -29, -29, -29,
            -29, -29, -29, -29, -29, -29, -29, -29, -28, -28, -28, -28, -28,
            -28, -23, -23, -23, -23, -23, -23, -23, -23, -23, -23, -23, -15,
            -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -14, -14, -14,
            -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14,
            -14, -14, -14, -14, -14, -14, -14, -14, -12, -12, -12, -12, -12,
            -12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -10, -10,
            -10, -10, -10, -10, -10, -10, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
            13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14,
            15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 18, 18, 18, 18, 18,
            18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
            18, 18, 18, 18, 18, 18, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
            20, 20, 20, 20, 20, 20, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
            25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 29, 29, 29, 29, 29, 29,
            29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
            30, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36,
            36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36,
            39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
            45, 45, 45, 45, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 48, 48,
            54, 54, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, 55, 55, 55, 55,
            55, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 62, 62, 62,
            62, 62, 62, 64, 64, 64, 64, 64, 64, 64, 64, 64, 66, 66, 66, 66, 66,
            66, 66, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 70,
            70, 70, 70, 70, 70, 72, 72, 72, 72, 72, 72, 72, 74, 74, 74, 74, 74,
            74, 74, 74, 74, 76, 76, 76, 76, 76, 76, 77, 77, 77, 77, 77, 77, 77,
            79, 79, 79, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 83, 83, 83, 83,
            83, 83, 83, 83, 83, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 98, 98,
            98, 98, 98, 98, 98, 98, 98, 98, 98, 102, 102, 102, 105, 105, 105,
            105, 105, 105, 105, 108, 108, 108, 108, 108, 108, 108, 108, 108,
            108, 114, 114, 114, 114, 114, 114, 114, 114, 120, 120, 120, 120,
            120, 120, 120, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121,
            121, 121, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126,
            126, 126, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 138,
            138, 138, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139,
            139, 139, 139, 148, 148, 148, 148, 148, 148, 152, 152, 152, 152,
            152, 152, 152, 152, 152, 152, 152, 164, 164, 164, 164, 164, 164,
            164, 164, 164, 164, 164, 164, 164, 164, 165, 165, 165, 165, 165,
            165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 167, 167,
            167, 167, 167, 167, 167, 167, 167, 172, 172, 172, 172, 172, 172,
            172, 172, 172, 172, 172, 172, 172, 172, 172, 172, 172, 172, 172,
            172, 172, 172, 172, 172, 172, 172, 172, 172, 172, 176, 176, 176,
            176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 178, 178, 178,
            178, 178, 178, 178, 178, 178, 178, 178, 178, 184, 184, 184, 184,
            184, 184, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185,
            185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185,
            185, 185, 185, 192, 192, 192, 192, 192, 192, 192, 192, 192, 203,
            203, 203, 203, 203, 207, 207, 207, 207, 207, 207, 207, 207, 207,
            207, 211, 211, 211, 211, 211, 211, 211, 211, 216, 216, 216, 216,
            216, 216, 216, 223, 223, 223, 223, 223, 223, 223, 223, 223, 223,
            223, 228, 228, 228, 228, 228, 228, 228, 228, 228, 229, 229, 229,
            229, 229, 229, 229, 229, 229, 229, 229, 229, 233, 233, 233, 233,
            233, 233, 241, 241, 241, 241, 241, 241, 241, 241, 241, 241, 241,
            241, 241, 247, 247, 247, 247, 247, 247, 247, 256, 256, 256, 256,
            256, 256, 256, 262, 262, 262, 262, 262]

        self._checkPositiveNegative(295, values)

    def testBinCountsPostiveNegative2(self):
        """
        See docstring for self._checkPositiveNegative (above).
        """
        self._checkPositiveNegative(11, [-262, 0, 100])

    def testBinCountsPostiveNegative3(self):
        """
        See docstring for self._checkPositiveNegative (above).
        """
        self._checkPositiveNegative(7, [-10, 0, 1, 1, 10])

    def testBinCountsPostiveNegative4(self):
        """
        See docstring for self._checkPositiveNegative (above).
        """
        self._checkPositiveNegative(27, [-5, 0, 3, 3, 5])

    def testGetItem(self):
        """
        The __getitem__ method must return the correct bin.
        """
        h = Histogram(3)
        list(map(h.add, range(9)))
        h.finalize()
        self.assertEqual([0, 1, 2], h[0])
        self.assertEqual([3, 4, 5], h[1])

    def testGetItemInvalidIndex(self):
        """
        The __getitem__ method must raise IndexError if passed the index
        of a non-existent bin.
        """
        h = Histogram(3)
        list(map(h.add, range(9)))
        h.finalize()
        six.assertRaisesRegex(self, IndexError, '^list index out of range$',
                              h.__getitem__, 4)

    def testHeightNoValues(self):
        """
        If a histogram with no value is finalized the height of each bin must
        be 0.
        """
        h = Histogram(3)
        h.finalize()
        for bin_ in h.bins:
            self.assertEqual(0, bin_.height)

    def testHeightWithValuesNoHeightSpecified(self):
        """
        If a histogram with added values but without specified heights is
        finalized, the height for each bin must be equal to the number of
        values in each bin.
        """
        h = Histogram(3)
        h.add(0)
        h.add(1)
        h.add(1)
        h.add(1)
        h.finalize()
        self.assertEqual(h.bins[0].height, len(h.bins[0]))
        self.assertEqual(h.bins[1].height, len(h.bins[1]))
        self.assertEqual(h.bins[2].height, len(h.bins[2]))

    def testHeightWithValuesHeightSpecified(self):
        """
        If a histogram with added values and with specified heights is
        finalized, the correct height for each bin must be calculated.
        """
        h = Histogram(3)
        h.add(0, height=3)
        h.add(1, height=1)
        h.add(1, height=2)
        h.add(1, height=3)
        h.add(2)
        h.finalize()
        self.assertEqual(3, h.bins[0].height)
        self.assertEqual(6, h.bins[1].height)
        self.assertEqual(1, h.bins[2].height)
