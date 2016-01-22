from unittest import TestCase

from light.histogram import Histogram
from light.overall_score import OverallScore


class TestOverallScore(TestCase):
    """
    Tests for the light.overall_score.OverallScore class.
    """
    def testEmptyHistogram(self):
        """
        An empty histogram must return a score of 0.0.
        """
        histogram = Histogram()
        histogram.finalize()
        os = OverallScore(histogram, [])
        score, analysis = os.calculateScore()
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'score': score,
            },
            analysis)

    def testPrintAnalysis(self):
        """
        The analysis of the overall score calculation must print correctly.
        """
        histogram = Histogram()
        histogram.finalize()
        os = OverallScore(histogram, [])
        score, analysis = os.calculateScore()

        expected = (
            'Overall score\nScore: %.4f' % score)

        self.assertEqual(expected, OverallScore.printAnalysis(analysis))
