from unittest import TestCase

from light.features import Landmark, TrigPoint
from light.histogram import Histogram
from light.significance import (Always, HashFraction, MaxBinHeight,
                                MeanBinHeight, AAFraction, getHeight)

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


class TestGetHeight(TestCase):
    """
    Tests for the light.significance.getHeight function.
    """
    def testGetHeightEmptyBin(self):
        """
        If an empty bin is passed, 0 must be returned.
        """
        length = getHeight([])
        self.assertEqual(0, length)

    def testGetHeightOnePair(self):
        """
        If the bin only contains one pair, the correct height must be returned.
        """
        bin_ = [{
            'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
            'subjectTrigPoint': TrigPoint('Peaks', 'P', 21),
            'queryLandmark': Landmark('AlphaHelix', 'A', 10, 9),
            'queryTrigPoint': TrigPoint('Peaks', 'P', 25),
        }]
        length = getHeight(bin_)
        self.assertEqual(20, length)

    def testGetHeightNoOverlap(self):
        """
        If there are no overlapping features in the subject and the query, the
        correct height must be returned.
        """
        bin_ = [{
            'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
            'subjectTrigPoint': TrigPoint('Peaks', 'P', 21),
            'queryLandmark': Landmark('AlphaHelix', 'A', 10, 9),
            'queryTrigPoint': TrigPoint('Peaks', 'P', 25),
        }, {
            'subjectLandmark': Landmark('AlphaHelix', 'A', 30, 9),
            'subjectTrigPoint': TrigPoint('Peaks', 'P', 41),
            'queryLandmark': Landmark('AlphaHelix', 'A', 40, 9),
            'queryTrigPoint': TrigPoint('Peaks', 'P', 55),
        }]
        length = getHeight(bin_)
        self.assertEqual(40, length)

    def testGetHeightOverlap(self):
        """
        If there is overlap between the features in the query, the correct
        height must be returned.
        """
        bin_ = [{
            'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
            'subjectTrigPoint': TrigPoint('Peaks', 'P', 21),
            'queryLandmark': Landmark('AlphaHelix', 'A', 10, 9),
            'queryTrigPoint': TrigPoint('Peaks', 'P', 25),
        }, {
            'subjectLandmark': Landmark('AlphaHelix', 'A', 5, 9),
            'subjectTrigPoint': TrigPoint('Peaks', 'P', 26),
            'queryLandmark': Landmark('AlphaHelix', 'A', 15, 9),
            'queryTrigPoint': TrigPoint('Peaks', 'P', 30),
        }]
        length = getHeight(bin_)
        self.assertEqual(32, length)


class TestAAFraction(TestCase):
    """
    Tests for the light.significance.AAFraction class.
    """
    def testAAFractionWhenNotSignificant(self):
        """
        The isSignificant method must return False if asked about a bin that is
        not significant.
        """
        match = {
            'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
            'subjectTrigPoint': TrigPoint('Peaks', 'P', 21),
            'queryLandmark': Landmark('AlphaHelix', 'A', 10, 9),
            'queryTrigPoint': TrigPoint('Peaks', 'P', 25),
        }
        histogram = Histogram(3)
        histogram.add(0, match)
        histogram.add(1, match)
        histogram.add(2, match)
        histogram.finalize()
        significance = AAFraction(histogram, 100, 0.75)
        self.assertFalse(significance.isSignificant(0))

    def testAAFractionWhenSignificant(self):
        """
        The isSignificant method must return True if asked about a bin that is
        significant.
        """
        match = {
            'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
            'subjectTrigPoint': TrigPoint('Peaks', 'P', 21),
            'queryLandmark': Landmark('AlphaHelix', 'A', 10, 9),
            'queryTrigPoint': TrigPoint('Peaks', 'P', 25),
        }
        histogram = Histogram(3)
        histogram.add(0, match)
        histogram.add(1, match)
        histogram.add(2, match)
        histogram.finalize()
        print(histogram.bins)
        significance = AAFraction(histogram, 10, 0.75)
        self.assertTrue(significance.isSignificant(0))

    def testAAFractionSignificanceAnalysis(self):
        """
        The correct analysis must be provided.
        """
        match = {
            'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
            'subjectTrigPoint': TrigPoint('Peaks', 'P', 21),
            'queryLandmark': Landmark('AlphaHelix', 'A', 10, 9),
            'queryTrigPoint': TrigPoint('Peaks', 'P', 25),
        }
        histogram = Histogram(3)
        histogram.add(0, match)
        histogram.add(1, match)
        histogram.add(2, match)
        histogram.finalize()
        significance = AAFraction(histogram, 10, 0.75)
        self.assertTrue(significance.isSignificant(0))
        self.assertEqual({'significanceCutoff': 7.5,
                          'significanceMethod': 'AAFraction'},
                         significance.analysis)
