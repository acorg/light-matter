from unittest import TestCase

from light.histogram import Histogram
from light.significance import Significance


class TestSignificance(TestCase):
    """
    Tests for the light.significance.Significance class.
    """
    def testHashFractionSignificantBins(self):
        """
        If a histogram with significant bins is passed, significant bins must
        be returned.
        """
        bins = [1, 1, 1, 1, 1, 6, 7, 8, 9]
        histogram = Histogram(5)
        for bin_ in bins:
            histogram.add(bin_)
        histogram.finalize()
        significanceInfo = Significance(histogram)
        signBins = significanceInfo.hashFraction(0.5, 5)
        self.assertEqual([{'bin': [1, 1, 1, 1, 1], 'index': 0, 'score': 1.0}],
                         signBins)

    def testHashFractionNoSignificantBins(self):
        """
        If a histogram without significant bins is passed, no significant bins
        must be returned.
        """
        bins = [1, 2, 3, 4, 5]
        histogram = Histogram(5)
        for bin_ in bins:
            histogram.add(bin_)
        histogram.finalize()
        significanceInfo = Significance(histogram)
        signBins = significanceInfo.hashFraction(0.3, 5)
        self.assertEqual([], signBins)

    def testBestScoreSignificantBins(self):
        """
        If no bins are significant, the best score must be None.
        """
        bins = [1, 1, 1, 1, 1, 5, 6, 8, 9, 10]
        histogram = Histogram(5)
        for bin_ in bins:
            histogram.add(bin_)
        histogram.finalize()
        significanceInfo = Significance(histogram)
        significanceInfo.hashFraction(0.1, 5)
        bestScore = significanceInfo.getBestScore()
        self.assertEqual(1.0, bestScore)

    def testBestScoreNoSignificantBins(self):
        """
        If there are significant bins, the right best score must be returned.
        """
        bins = [1, 2, 3, 4, 5]
        histogram = Histogram(5)
        for bin_ in bins:
            histogram.add(bin_)
        histogram.finalize()
        significanceInfo = Significance(histogram)
        significanceInfo.hashFraction(0.3, 5)
        bestScore = significanceInfo.getBestScore()
        self.assertEqual(None, bestScore)
