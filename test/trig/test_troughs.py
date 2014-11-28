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
        troughs = Troughs()
        result = troughs.convertAAToProperties('ASDGEAHSDTDSCV')
        self.assertEqual([-1.401365904912, -0.17713381906999998,
                         0.003636363636359996, -0.6484712512417,
                         -0.391107796081, -1.401365904912,
                         -0.013648991654999942, -0.17713381906999998,
                         0.003636363636359996, -0.7214954158464001,
                         0.003636363636359996, -0.17713381906999998,
                         -0.2761322022899999, -1.9545882971], result)

    def testFindWithoutTrough(self):
        """
        The find method must return an empty generator when no trough is
        present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        troughs = Troughs()
        result = list(troughs.find(read))
        self.assertEqual([], result)

    def testTwoTroughs(self):
        """
        The find method must find two troughs.
        """
        read = AARead('id', 'AVAAAAAVAAA')
        troughs = Troughs()
        result = list(troughs.find(read))
        self.assertEqual([TrigPoint('Troughs', 'T', 1),
                          TrigPoint('Troughs', 'T', 7)], result)

    def testNoTroughAtBeginning(self):
        """
        The find method must not find a troughs at the beginning of the
        sequence.
        """
        read = AARead('id', 'VAAAAAAAAA')
        troughs = Troughs()
        result = list(troughs.find(read))
        self.assertEqual([], result)

    def testNoTroughAtEnd(self):
        """
        The find method must not find a trough at the end of the sequence.
        """
        read = AARead('id', 'AAAAAAV')
        troughs = Troughs()
        result = list(troughs.find(read))
        self.assertEqual([], result)
