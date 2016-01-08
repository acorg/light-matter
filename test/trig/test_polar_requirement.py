from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.polar_requirement import PolarRequirement


class TestPolarRequirement(TestCase):
    """
    Tests for the trig.polar_requirement.PolarRequirement class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right cluster list.
        """
        pr = PolarRequirement()
        result = pr.convertAAToPolarRequirementClusters('ASDGEAHSDTDSCV')
        self.assertEqual([2, 2, 4, 2, 4, 2, 2, 2, 4, 2, 4, 2, 1, 1], result)

    def testConversionAANotInProperties(self):
        """
        If an aa is not present in the dark.aa.PROPERTY_CLUSTERS dict, the
        sequence must be converted to the right cluster list.
        """
        pr = PolarRequirement()
        result = pr.convertAAToPolarRequirementClusters('ASDGEAHSDTDSCVX')
        self.assertEqual([2, 2, 4, 2, 4, 2, 2, 2, 4, 2, 4, 2, 1, 1, 0], result)

    def testFindWithoutPolarRequirementPeak(self):
        """
        The find method must return an empty generator when no polar
        requirement peak is present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        polarRequirement = PolarRequirement()
        result = list(polarRequirement.find(read))
        self.assertEqual([], result)

    def testTwoPolarRequirementPeaks(self):
        """
        The find method must find two polar requirement peaks.
        """
        read = AARead('id', 'GIGGGGGIGGG')
        polarRequirement = PolarRequirement()
        result = list(polarRequirement.find(read))
        self.assertEqual([TrigPoint('PolarRequirement', 'PR', 1),
                          TrigPoint('PolarRequirement', 'PR', 7)], result)

    def testNoPolarRequirementPeakAtBeginning(self):
        """
        The find method must not find a polar requirement peak at the beginning
        of the sequence.
        """
        read = AARead('id', 'GIIIIIII')
        polarRequirement = PolarRequirement()
        result = list(polarRequirement.find(read))
        self.assertEqual([], result)

    def testNoPolarRequirementPeakAtEnd(self):
        """
        The find method must not find a polar requirement peak at the end of
        the sequence.
        """
        read = AARead('id', 'GGGGGGI')
        polarRequirement = PolarRequirement()
        result = list(polarRequirement.find(read))
        self.assertEqual([], result)

    def testMixedSequence(self):
        """
        The right polar requirement peaks must be found in a mixed sequence.
        """
        read = AARead('id', 'GHGASFQLLDFGF')
        polarRequirement = PolarRequirement()
        result = list(polarRequirement.find(read))
        self.assertEqual([TrigPoint('PolarRequirement', 'PR', 5),
                          TrigPoint('PolarRequirement', 'PR', 6),
                          TrigPoint('PolarRequirement', 'PR', 9),
                          TrigPoint('PolarRequirement', 'PR', 11)], result)
