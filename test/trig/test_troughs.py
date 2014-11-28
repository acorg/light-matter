from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.troughs import Troughs


class TestTroughs(TestCase):
    """
    Tests for the TrigPoint.Troughs class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right properties
        string.
        """
        trigPoint = Troughs()
        result = trigPoint.convertAAToProperties('ASDGEAHSDTDSCV')
        self.assertEqual([-1.401365904912, -0.17713381906999998,
                         0.003636363636359996, -0.6484712512417,
                         -0.391107796081, -1.401365904912,
                         -0.013648991654999942, -0.17713381906999998,
                         0.003636363636359996, -0.7214954158464001,
                         0.003636363636359996, -0.17713381906999998,
                         -0.2761322022899999, -1.9545882971], result)

    def testFindWithoutPeak(self):
        """
        The find method must return an empty generator when no peak is
        present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        trigPoint = Troughs()
        result = list(trigPoint.find(read))
        self.assertEqual([], result)

    def testTwoPeaks(self):
        """
        The find method must find two peaks.
        """
        read = AARead('id', 'AVAAAAAVAAA')
        trigPoint = Troughs()
        result = list(trigPoint.find(read))
        self.assertEqual([TrigPoint('Troughs', 'T', 1),
                          TrigPoint('Troughs', 'T', 7)], result)

    def testNoPeakAtBeginning(self):
        """
        The find method must not find a peak at the beginning of the sequence.
        """
        read = AARead('id', 'VAAAAAAAAA')
        trigPoint = Troughs()
        result = list(trigPoint.find(read))
        self.assertEqual([], result)

    def testNoPeakAtEnd(self):
        """
        The find method must not find a peak at the end of the sequence.
        """
        read = AARead('id', 'AAAAAAV')
        trigPoint = Troughs()
        result = list(trigPoint.find(read))
        self.assertEqual([], result)
