from unittest import TestCase

from light.histogram import Histogram
from light.significance import HashFraction


class TestHashFraction(TestCase):
    """
    Tests for the light.significance.HashFraction class.
    """
    def testIsSignificantWhenNotSignificant(self):
        """
        The isSignificant method must return False if asked about a bin
        that is not significant.
        """
        histogram = Histogram(5)
        map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9])
        histogram.finalize()
        significance = HashFraction(histogram, 10, 0.75)
        self.assertFalse(significance.isSignificant(0))

    def testIsSignificantWhenSignificant(self):
        """
        The isSignificant method must return True if asked about a bin
        that is significant.
        """
        histogram = Histogram(5)
        map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9])
        histogram.finalize()
        significance = HashFraction(histogram, 10, 0.1)
        self.assertTrue(significance.isSignificant(0))
