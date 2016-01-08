from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.hydropathy import Hydropathy


class TestHydropathy(TestCase):
    """
    Tests for the trig.Hydropathy class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right cluster
        list.
        """
        hydropathy = Hydropathy()
        result = hydropathy.convertAAToHydropathyClusters('ASDGEAHSDTDSCV')
        self.assertEqual([3, 2, 1, 2, 1, 3, 1, 2, 1, 2, 1, 2, 3, 4], result)

    def testConversionAANotInProperties(self):
        """
        If an aa is not present in the dark.aa.PROPERTY_CLUSTERS dict, the
        sequence must be converted to the right cluster list.
        """
        hydropathy = Hydropathy()
        result = hydropathy.convertAAToHydropathyClusters('ASDGEAHSDTDSCVX')
        self.assertEqual([3, 2, 1, 2, 1, 3, 1, 2, 1, 2, 1, 2, 3, 4, 0], result)

    def testFindWithoutHydropathyPeak(self):
        """
        The find method must return an empty generator when no hydropathy peak
        is present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        hydropathy = Hydropathy()
        result = list(hydropathy.find(read))
        self.assertEqual([], result)

    def testTwoHydropathyPeaks(self):
        """
        The find method must find two hydropathy peaks.
        """
        read = AARead('id', 'GHGGGGGHGGG')
        hydropathy = Hydropathy()
        result = list(hydropathy.find(read))
        self.assertEqual([TrigPoint('Hydropathy', 'H', 1),
                          TrigPoint('Hydropathy', 'H', 7)], result)

    def testNoHydropathyPeakAtBeginning(self):
        """
        The find method must not find a hydropathy peak at the beginning of the
        sequence.
        """
        read = AARead('id', 'GHHHHHHHH')
        hydropathy = Hydropathy()
        result = list(hydropathy.find(read))
        self.assertEqual([], result)

    def testNoHydropathyPeakAtEnd(self):
        """
        The find method must not find a Hydropathy peak at the end of the
        sequence.
        """
        read = AARead('id', 'GGGGGGH')
        hydropathy = Hydropathy()
        result = list(hydropathy.find(read))
        self.assertEqual([], result)

    def testMixedSequence(self):
        """
        The right hydropathy peaks must be found in a mixed sequence.
        """
        read = AARead('id', 'GHGASFQLLDFGF')
        hydropathy = Hydropathy()
        result = list(hydropathy.find(read))
        self.assertEqual([TrigPoint('Hydropathy', 'H', 1),
                          TrigPoint('Hydropathy', 'H', 3),
                          TrigPoint('Hydropathy', 'H', 4),
                          TrigPoint('Hydropathy', 'H', 11)], result)
