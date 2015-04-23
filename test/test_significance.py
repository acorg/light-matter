from unittest import TestCase

from light.histogram import Histogram
from light.significance import Always, HashFraction, MaxBinHeight

from test.sample_data import DB, COWPOX


class TestHashFraction(TestCase):
    """
    Tests for the light.significance.HashFraction class.
    """
    def testHashFractionIsSignificantWhenNotSignificant(self):
        """
        The isSignificant method must return False if asked about a bin
        that is not significant.
        """
        histogram = Histogram(5)
        map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9])
        histogram.finalize()
        significance = HashFraction(histogram, 10, 0.75)
        self.assertFalse(significance.isSignificant(0))

    def testHashFractionIsSignificantWhenSignificant(self):
        """
        The isSignificant method must return True if asked about a bin
        that is significant.
        """
        histogram = Histogram(5)
        map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9])
        histogram.finalize()
        significance = HashFraction(histogram, 10, 0.1)
        self.assertTrue(significance.isSignificant(0))


class TestMaxBinHeight(TestCase):
    """
    Tests for the light.significance.MaxBinHeight class.
    """
    def testMaxBinHeightIsSignificantWhenNotSignificant(self):
        """
        The isSignificant method must return False if asked about a bin
        that is not significant.
        """
        histogram = Histogram(5)
        map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9])
        histogram.finalize()
        significance = MaxBinHeight(histogram, COWPOX, DB)
        self.assertFalse(significance.isSignificant(1))

    def testMaxBinHeightIsSignificantWhenSignificant(self):
        """
        The isSignificant method must return True if asked about a bin
        that is significant.
        """
        histogram = Histogram(5)
        map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9])
        histogram.finalize()
        significance = MaxBinHeight(histogram, COWPOX, DB)
        self.assertTrue(significance.isSignificant(0))


class TestAlways(TestCase):
    """
    Tests for the light.significance.Always class.
    """
    def testAlwaysMustBeTrue(self):
        """
        The Always significance method must return True.
        """
        histogram = Histogram(5)
        map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9])
        histogram.finalize()
        significance = Always()
        self.assertTrue(significance.isSignificant(0))
