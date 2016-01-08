from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.composition import Composition


class TestComposition(TestCase):
    """
    Tests for the trig.Composition class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right cluster
        list.
        """
        composition = Composition()
        result = composition.convertAAToCompositionClusters('ASDGEAHSDTDSCV')
        self.assertEqual([1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 2, 3, 1], result)

    def testConversionAANotInProperties(self):
        """
        If an aa is not present in the dark.aa.PROPERTY_CLUSTERS dict, the
        sequence must be converted to the right cluster list.
        """
        composition = Composition()
        result = composition.convertAAToCompositionClusters('ASDGEAHSDTDSCVX')
        self.assertEqual([1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 2, 3, 1, 0], result)

    def testFindWithoutCompositionPeak(self):
        """
        The find method must return an empty generator when no composition peak
        is present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        composition = Composition()
        result = list(composition.find(read))
        self.assertEqual([], result)

    def testTwoCompositionPeaks(self):
        """
        The find method must find two composition peaks.
        """
        read = AARead('id', 'GNGGGGGCGGG')
        composition = Composition()
        result = list(composition.find(read))
        self.assertEqual([TrigPoint('Composition', 'CP', 1),
                          TrigPoint('Composition', 'CP', 7)], result)

    def testNoCompositionPeakAtBeginning(self):
        """
        The find method must not find a composition peak at the beginning of
        the sequence.
        """
        read = AARead('id', 'GHHHHHHHH')
        composition = Composition()
        result = list(composition.find(read))
        self.assertEqual([], result)

    def testNoCompositionPeakAtEnd(self):
        """
        The find method must not find a composition peak at the end of the
        sequence.
        """
        read = AARead('id', 'GGGGGGH')
        composition = Composition()
        result = list(composition.find(read))
        self.assertEqual([], result)

    def testMixedSequence(self):
        """
        The right composition peaks must be found in a mixed sequence.
        """
        read = AARead('id', 'GHGASFQLLDFGF')
        composition = Composition()
        result = list(composition.find(read))
        self.assertEqual([TrigPoint('Composition', 'CP', 4),
                          TrigPoint('Composition', 'CP', 9)], result)
