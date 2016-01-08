from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.aromaticity import Aromaticity


class TestAromaticity(TestCase):
    """
    Tests for the trig.Aromaticity class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right cluster
        list.
        """
        aromaticity = Aromaticity()
        result = aromaticity.convertAAToAromaticityClusters('ASDWEAHYDTDSCV')
        self.assertEqual([1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1], result)

    def testConversionAANotInProperties(self):
        """
        If an aa is not present in the dark.aa.PROPERTY_CLUSTERS dict, the
        sequence must be converted to the right cluster list.
        """
        aromaticity = Aromaticity()
        result = aromaticity.convertAAToAromaticityClusters('ASDWEAHYDTDSCVX')
        self.assertEqual([1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 0], result)

    def testFindWithoutAromaticityPeak(self):
        """
        The find method must return an empty generator when no aromaticity peak
        is present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        aromaticity = Aromaticity()
        result = list(aromaticity.find(read))
        self.assertEqual([], result)

    def testTwoAromaticityPeaks(self):
        """
        The find method must find two aromaticity peaks.
        """
        read = AARead('id', 'GYGGGGGYGGG')
        aromaticity = Aromaticity()
        result = list(aromaticity.find(read))
        self.assertEqual([TrigPoint('Aromaticity', 'AR', 1),
                          TrigPoint('Aromaticity', 'AR', 7)], result)

    def testNoAromaticityPeakAtBeginning(self):
        """
        The find method must not find a aromaticity peak at the beginning of
        the sequence.
        """
        read = AARead('id', 'GHHHHHHHH')
        aromaticity = Aromaticity()
        result = list(aromaticity.find(read))
        self.assertEqual([], result)

    def testNoAromaticityPeakAtEnd(self):
        """
        The find method must not find a aromaticity peak at the end of the
        sequence.
        """
        read = AARead('id', 'GGGGGGH')
        aromaticity = Aromaticity()
        result = list(aromaticity.find(read))
        self.assertEqual([], result)

    def testMixedSequence(self):
        """
        The right aromaticity peaks must be found in a mixed sequence.
        """
        read = AARead('id', 'GHGASFQLLDFGF')
        aromaticity = Aromaticity()
        result = list(aromaticity.find(read))
        self.assertEqual([TrigPoint('Aromaticity', 'AR', 1),
                          TrigPoint('Aromaticity', 'AR', 5),
                          TrigPoint('Aromaticity', 'AR', 10),
                          TrigPoint('Aromaticity', 'AR', 11)], result)
