from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.hydroxythiolation import Hydroxythiolation


class TestHydroxythiolation(TestCase):
    """
    Tests for the trig.Hydroxythiolation class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right cluster list.
        """
        hyd = Hydroxythiolation()
        result = hyd.convertAAToHydroxythiolationClusters('ASDGEAHSDTDSCV')
        self.assertEqual([2, 5, 2, 2, 1, 2, 3, 5, 2, 5, 2, 5, 5, 1], result)

    def testConversionAANotInProperties(self):
        """
        If an aa is not present in the dark.aa.PROPERTY_CLUSTERS dict, the
        sequence must be converted to the right cluster list.
        """
        hyd = Hydroxythiolation()
        result = hyd.convertAAToHydroxythiolationClusters('ASDGEAHSDTDSCVX')
        self.assertEqual([2, 5, 2, 2, 1, 2, 3, 5, 2, 5, 2, 5, 5, 1, 0], result)

    def testFindWithoutHydroxythiolationPeak(self):
        """
        The find method must return an empty generator when no
        hydroxythiolation peak is present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        hyd = Hydroxythiolation()
        result = list(hyd.find(read))
        self.assertEqual([], result)

    def testTwoHydroxythiolationPeaks(self):
        """
        The find method must find two hydroxythiolation peaks.
        """
        read = AARead('id', 'GHGGGGGHGGG')
        hyd = Hydroxythiolation()
        result = list(hyd.find(read))
        self.assertEqual([TrigPoint('Hydroxythiolation', 'HT', 1),
                          TrigPoint('Hydroxythiolation', 'HT', 7)], result)

    def testNoHydroxythiolationPeakAtBeginning(self):
        """
        The find method must not find a hydroxythiolation peak at the beginning
        of the sequence.
        """
        read = AARead('id', 'GHHHHHHHH')
        hyd = Hydroxythiolation()
        result = list(hyd.find(read))
        self.assertEqual([], result)

    def testNoHydroxythiolationPeakAtEnd(self):
        """
        The find method must not find a hydroxythiolation peak at the end of
        the sequence.
        """
        read = AARead('id', 'GGGGGGH')
        hyd = Hydroxythiolation()
        result = list(hyd.find(read))
        self.assertEqual([], result)

    def testMixedSequence(self):
        """
        The right hydroxythiolation peaks must be found in a mixed sequence.
        """
        read = AARead('id', 'GHGASFQLLDFGF')
        hyd = Hydroxythiolation()
        result = list(hyd.find(read))
        self.assertEqual([TrigPoint('Hydroxythiolation', 'HT', 1),
                          TrigPoint('Hydroxythiolation', 'HT', 10),
                          TrigPoint('Hydroxythiolation', 'HT', 11)], result)
