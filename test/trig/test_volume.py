from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.volume import Volume


class TestVolume(TestCase):
    """
    Tests for the TrigPoint.Volume class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right cluster
        list.
        """
        volume = Volume()
        result = volume.convertAAToVolumeClusters('ASDGEAHSDTDSCV')
        self.assertEqual([2, 2, 3, 1, 4, 2, 4, 2, 3, 3, 3, 2, 3, 4], result)

    def testConversionAANotInProperties(self):
        """
        If an aa is not present in the dark.aa.PROPERTY_CLUSTERS dict, the
        sequence must be converted to the right cluster list.
        """
        volume = Volume()
        result = volume.convertAAToVolumeClusters('ASDGEAHSDTDSCVX')
        self.assertEqual([2, 2, 3, 1, 4, 2, 4, 2, 3, 3, 3, 2, 3, 4, 0], result)

    def testFindWithoutVolumePeak(self):
        """
        The find method must return an empty generator when no volume peak is
        present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        volume = Volume()
        result = list(volume.find(read))
        self.assertEqual([], result)

    def testTwoVolumePeaks(self):
        """
        The find method must find two volume peaks.
        """
        read = AARead('id', 'GHGGGGGHGGG')
        volume = Volume()
        result = list(volume.find(read))
        self.assertEqual([TrigPoint('Volume', 'V', 1),
                          TrigPoint('Volume', 'V', 7)], result)

    def testNoVolumePeakAtBeginning(self):
        """
        The find method must not find a volume peak at the beginning of the
        sequence.
        """
        read = AARead('id', 'GHHHHHHHH')
        volume = Volume()
        result = list(volume.find(read))
        self.assertEqual([], result)

    def testNoVolumePeakAtEnd(self):
        """
        The find method must not find a volume peak at the end of the sequence.
        """
        read = AARead('id', 'GGGGGGH')
        volume = Volume()
        result = list(volume.find(read))
        self.assertEqual([], result)

    def testMixedSequence(self):
        """
        The right volume peaks must be found in a mixed sequence.
        """
        read = AARead('id', 'GHGASFQLLDFGF')
        volume = Volume()
        result = list(volume.find(read))
        self.assertEqual([TrigPoint('Volume', 'V', 1),
                          TrigPoint('Volume', 'V', 9),
                          TrigPoint('Volume', 'V', 11)], result)
