import warnings
from unittest import TestCase

from light.score import MinHashesScore, FeatureMatchingScore
from light.histogram import Histogram


class TestMinHashesScore(TestCase):
    """
    Tests for the light.score.MinHashesScore class.
    """
    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of 0.0
        """
        histogram = Histogram()
        histogram.finalize()
        mhs = MinHashesScore(histogram, 0)
        self.assertEqual(0.0, mhs.calculateScore(0))

    def testZeroMinHashCount(self):
        """
        If the minimum hash count for the query/subject is zero, the score
        for any bin must be 0.0
        """
        # Use a histogram with one bin, and put something in it, to make sure
        # we're not testing an empty histogram bin.
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 0)
        self.assertEqual(0.0, mhs.calculateScore(0))

    def testBinWithOneElementMinHashCountOne(self):
        """
        If the minimum hash count for the query/subject is one and the bin has
        one thing in it, the score for the bin must be 1.0
        """
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 1)
        self.assertEqual(1.0, mhs.calculateScore(0))

    def testBinWithOneElementMinHashCountTwo(self):
        """
        If the minimum hash count for the query/subject is two and the bin has
        one thing in it, the score for the bin must be 0.5
        """
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 2)
        self.assertEqual(0.5, mhs.calculateScore(0))

    def testBinWithTwoElementsMinHashCountOne(self):
        """
        If the minimum hash count for the query/subject is one and the bin has
        two things in it, a warning will be issued to say that the bin has more
        things in it than should be possible (given the minimum hash count of
        the query/subject).
        """
        histogram = Histogram(1)
        histogram.add(444)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 1)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            mhs.calculateScore(0)
            self.assertEqual(1, len(w))
            self.assertTrue(issubclass(w[0].category, RuntimeWarning))
            error = ('Bin contains 2 deltas for a query/subject pair with a '
                     'minimum hash count of only 1.')
            self.assertIn(error, str(w[0].message))


class TestFeatureMatchingScore(TestCase):
    """
    Tests for the light.score.FeatureMatchingScore class.
    """
    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of 0.0
        """
        histogram = Histogram()
        histogram.finalize()
        mhs = FeatureMatchingScore(histogram)
        self.assertEqual(0.0, mhs.calculateScore(0))
