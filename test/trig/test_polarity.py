from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.polarity import Polarity


class TestPolarity(TestCase):
    """
    Tests for the trig.polarity.Polarity class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right cluster list.
        """
        polarity = Polarity()
        result = polarity.convertAAToPolarityClusters('ASDGEAHSDTDSCV')
        self.assertEqual([2, 2, 4, 2, 4, 2, 3, 2, 4, 2, 4, 2, 1, 1], result)

    def testConversionAANotInProperties(self):
        """
        If an aa is not present in the dark.aa.PROPERTY_CLUSTERS dict, the
        sequence must be converted to the right cluster list.
        """
        polarity = Polarity()
        result = polarity.convertAAToPolarityClusters('ASDGEAHSDTDSCVX')
        self.assertEqual([2, 2, 4, 2, 4, 2, 3, 2, 4, 2, 4, 2, 1, 1, 0], result)

    def testFindWithoutPolarityPeak(self):
        """
        The find method must return an empty generator when no polarity peak is
        present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        polarity = Polarity()
        result = list(polarity.find(read))
        self.assertEqual([], result)

    def testTwoPolarityPeaks(self):
        """
        The find method must find two polarity peaks.
        """
        read = AARead('id', 'GHGGGGGHGGG')
        polarity = Polarity()
        result = list(polarity.find(read))
        self.assertEqual([TrigPoint('Polarity', 'PO', 1),
                          TrigPoint('Polarity', 'PO', 7)], result)

    def testNoPolarityPeakAtBeginning(self):
        """
        The find method must not find a polarity peak at the beginning of the
        sequence.
        """
        read = AARead('id', 'GHHHHHHHH')
        polarity = Polarity()
        result = list(polarity.find(read))
        self.assertEqual([], result)

    def testNoPolarityPeakAtEnd(self):
        """
        The find method must not find a polarity peak at the end of the
        sequence.
        """
        read = AARead('id', 'GGGGGGH')
        polarity = Polarity()
        result = list(polarity.find(read))
        self.assertEqual([], result)

    def testMixedSequence(self):
        """
        The right polarity peaks must be found in a mixed sequence.
        """
        read = AARead('id', 'GHGASFQLLDFGF')
        polarity = Polarity()
        result = list(polarity.find(read))
        self.assertEqual([TrigPoint('Polarity', 'PO', 1),
                         TrigPoint('Polarity', 'PO', 6),
                          TrigPoint('Polarity', 'PO', 9),
                          TrigPoint('Polarity', 'PO', 11)], result)
