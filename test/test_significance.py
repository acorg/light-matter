from unittest import TestCase

from light.histogram import Histogram
from light.significance import (Always, HashFraction, MaxBinHeight,
                                MeanBinHeight)

from test.sample_data import DB, COWPOX, SQUIRRELPOX


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
        list(map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9]))
        histogram.finalize()
        significance = HashFraction(histogram, 10, 0.75)
        self.assertFalse(significance.isSignificant(0))

    def testHashFractionIsSignificantWhenSignificant(self):
        """
        The isSignificant method must return True if asked about a bin
        that is significant.
        """
        histogram = Histogram(5)
        list(map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9]))
        histogram.finalize()
        significance = HashFraction(histogram, 10, 0.1)
        self.assertTrue(significance.isSignificant(0))

    def testHashFractionSignificanceAnalysis(self):
        """
        The correct analysis must be provided.
        """
        histogram = Histogram(5)
        list(map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9]))
        histogram.finalize()
        significance = HashFraction(histogram, 10, 0.1)
        self.assertEqual({'significanceCutoff': 1.0,
                          'significanceMethod': 'HashFraction'},
                         significance.analysis)


class TestMeanBinHeight(TestCase):
    """
    Tests for the light.significance.MeanBinHeight class.
    """
    def testMeanBinHeightIsSignificantWhenNotSignificant(self):
        """
        The isSignificant method must return False if asked about a bin
        that is not significant.
        """
        histogram = Histogram(5)
        list(map(histogram.add, [1, 1, 1, 1, 1, 7, 8, 9]))
        histogram.finalize()
        significance = MeanBinHeight(histogram, SQUIRRELPOX, DB)
        self.assertFalse(significance.isSignificant(1))

    def testMeanBinHeightIsSignificantWhenSignificant(self):
        """
        The isSignificant method must return True if asked about a bin
        that is significant.
        """
        histogram = Histogram(5)
        list(map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9]))
        histogram.finalize()
        significance = MeanBinHeight(histogram, COWPOX, DB)
        self.assertTrue(significance.isSignificant(0))

    def testMeanBinHeightSignificanceAnalysis(self):
        """
        The right analysis must be returned.
        """
        histogram = Histogram(5)
        list(map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9]))
        histogram.finalize()
        significance = MeanBinHeight(histogram, COWPOX, DB)
        self.assertEqual({'meanBinHeight': 0.0,
                          'significanceCutoff': 0.0,
                          'significanceMethod': 'MeanBinHeight',
                          'standardDeviation': 0.0},
                         significance.analysis)


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
        list(map(histogram.add, [1, 1, 1, 1, 1, 7, 8, 9]))
        histogram.finalize()
        significance = MaxBinHeight(histogram, SQUIRRELPOX, DB)
        self.assertFalse(significance.isSignificant(1))

    def testMaxBinHeightIsSignificantWhenSignificant(self):
        """
        The isSignificant method must return True if asked about a bin
        that is significant.
        """
        histogram = Histogram(5)
        list(map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9]))
        histogram.finalize()
        significance = MaxBinHeight(histogram, COWPOX, DB)
        self.assertTrue(significance.isSignificant(0))

    def testMaxBinHeightSignificanceAnalysis(self):
        """
        The correct analysis must be provided.
        """
        histogram = Histogram(5)
        list(map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9]))
        histogram.finalize()
        significance = MaxBinHeight(histogram, COWPOX, DB)
        self.assertEqual({'significanceCutoff': 0.0,
                          'significanceMethod': 'MaxBinHeight'},
                         significance.analysis)


class TestAlways(TestCase):
    """
    Tests for the light.significance.Always class.
    """
    def testAlwaysMustBeTrue(self):
        """
        The Always significance method must return True.
        """
        histogram = Histogram(5)
        list(map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9]))
        histogram.finalize()
        significance = Always()
        self.assertTrue(significance.isSignificant(0))

    def testAlwaysSignificanceAnalysis(self):
        """
        The correct analysis must be provided.
        """
        histogram = Histogram(5)
        list(map(histogram.add, [1, 1, 1, 1, 1, 6, 7, 8, 9]))
        histogram.finalize()
        significance = Always()
        self.assertEqual({'significanceMethod': 'Always'},
                         significance.analysis)
