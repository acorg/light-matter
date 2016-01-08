from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.isoelectric_point import IsoelectricPoint


class TestIsoelectricPoint(TestCase):
    """
    Tests for the trig.IsoelectricPoint class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right cluster list.
        """
        iep = IsoelectricPoint()
        result = iep.convertAAToIEPClusters('ASDGEAHSDTDSCV')
        self.assertEqual([2, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 2, 2], result)

    def testConversionAANotInProperties(self):
        """
        If an aa is not present in the dark.aa.PROPERTY_CLUSTERS dict, the
        sequence must be converted to the right cluster list.
        """
        iep = IsoelectricPoint()
        result = iep.convertAAToIEPClusters('ASDGEAHSDTDSCVX')
        self.assertEqual([2, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 2, 2, 0], result)

    def testFindWithoutIsoelectricPointPeak(self):
        """
        The find method must return an empty generator when no isoelectric
        point peak is present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        iep = IsoelectricPoint()
        result = list(iep.find(read))
        self.assertEqual([], result)

    def testTwoIsoelectricPointPeaks(self):
        """
        The find method must find two isoelectric point peaks.
        """
        read = AARead('id', 'GHGGGGGHGGG')
        iep = IsoelectricPoint()
        result = list(iep.find(read))
        self.assertEqual([TrigPoint('IsoelectricPoint', 'IEP', 1),
                          TrigPoint('IsoelectricPoint', 'IEP', 7)], result)

    def testNoIsoelectricPointPeakAtBeginning(self):
        """
        The find method must not find an isoelectric point peak at the
        beginning of the sequence.
        """
        read = AARead('id', 'GHHHHHHHH')
        iep = IsoelectricPoint()
        result = list(iep.find(read))
        self.assertEqual([], result)

    def testNoIsoelectricPointPeakAtEnd(self):
        """
        The find method must not find an isoelectric point peak at the end of
        the sequence.
        """
        read = AARead('id', 'GGGGGGH')
        iep = IsoelectricPoint()
        result = list(iep.find(read))
        self.assertEqual([], result)

    def testMixedSequence(self):
        """
        The right isoelectric point peaks must be found in a mixed sequence.
        """
        read = AARead('id', 'GHGASFQLLDFGF')
        iep = IsoelectricPoint()
        result = list(iep.find(read))
        self.assertEqual([TrigPoint('IsoelectricPoint', 'IEP', 1),
                          TrigPoint('IsoelectricPoint', 'IEP', 9)], result)
