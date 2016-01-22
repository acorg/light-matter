from unittest import TestCase

from light.histogram import Histogram
from light.overall_score import BestBinScore


class TestBestBinScore(TestCase):
    """
    Tests for the light.overall_score.BestBinScore class.
    """
    def testEmptyHistogram(self):
        """
        An empty histogram must return a score of 0.0.
        """
        histogram = Histogram()
        histogram.finalize()
        bestBinScore = BestBinScore(histogram, [])
        score, analysis = bestBinScore.calculateScore()
        self.assertEqual(None, score)
        self.assertEqual(
            {
                'score': score,
                'scoreClass': bestBinScore.__class__,
            },
            analysis)

    def testPrintAnalysis(self):
        """
        The analysis of the overall score calculation must print correctly.
        """
        histogram = Histogram()
        histogram.finalize()
        bestBinScore = BestBinScore(histogram, [])
        score, analysis = bestBinScore.calculateScore()

        expected = (
            'Overall score method: BestBinScore\nOverall score: %s' % score)

        self.assertEqual(expected, BestBinScore.printAnalysis(analysis))
