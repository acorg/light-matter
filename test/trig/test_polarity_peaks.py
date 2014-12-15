from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.polarity_peaks import PolarityPeaks


class TestPolarityPeaks(TestCase):
    """
    Tests for the TrigPoint.PolarityPeaks class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right properties
        string.
        """
        peaks = PolarityPeaks()
        result = peaks.sumProperties('ASDGEAHSDTDSCV')
        self.assertEqual([-0.20987654321, -0.1481481481483, 0.8518518518517,
                          0.864197530864, 1.691358024691, 1.481481481481,
                          1.839506172839, 1.9012345679007001,
                          2.9012345679007003, 2.8148148148143, 3.8148148148143,
                          3.876543209876, 3.024691358024, 2.271604938271],
                         result)

    def testFindPolarityPeaks(self):
        """
        The find method must find a peak.
        """
        read = AARead('i', 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDLLLLLLLLLL')
        peaks = PolarityPeaks()
        result = list(peaks.find(read, windowSize=43))
        self.assertEqual([TrigPoint('PolarityPeak', 'O', 33)], result)

    def testPolarityPeakAtEnd(self):
        """
        If all aa are the same, the polarity peak must be at the start of the
        sequence.
        """
        read = AARead('i', 'LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL')
        peaks = PolarityPeaks()
        result = list(peaks.find(read, windowSize=43))
        self.assertEqual(TrigPoint('PolarityPeak', 'O', 0), result[0])

    def testFindPolarityPeaksSmallerWindow(self):
        """
        The find method must find a peak, with a windowSize smaller than
        the sequence length.
        """
        read = AARead('i', 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDLLLLLLLLLL')
        peaks = PolarityPeaks()
        result = list(peaks.find(read, windowSize=41))
        self.assertEqual([TrigPoint('PolarityPeak', 'O', 33),
                          TrigPoint('PolarityPeak', 'O', 33),
                          TrigPoint('PolarityPeak', 'O', 33)], result)
