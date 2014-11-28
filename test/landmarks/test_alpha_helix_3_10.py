from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.alpha_helix_3_10 import AlphaHelix_3_10


class TestAlphaHelix_3_10(TestCase):
    """
    Tests for the Landmark.AlphaHelix_3_10 class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right properties
        string.
        """
        landmark = AlphaHelix_3_10()
        result = landmark.convertAAToHydrophobicHydrophilic('ASDGEAHSDTDSCV')
        self.assertEqual('OIIIIOOIIOIIOO', result)

    def testFindWithoutHelix(self):
        """
        The find method must return an empty generator when no helix is
        present.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        landmark = AlphaHelix_3_10()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testFindWithOneRepeat(self):
        """
        The find method must return an empty generator if there is only one
        instance of the OII pattern
        """
        read = AARead('id', 'FRR')
        landmark = AlphaHelix_3_10()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testHelixTwoRepeats(self):
        """
        The find method must find a helix with a repeat count of two.
        """
        read = AARead('id', 'FRRFRRF')
        landmark = AlphaHelix_3_10()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('AlphaHelix_3_10', 'B', 0, 7, 2)], result)

    def testHelixTwoRepeatsWithNonZeroOffset(self):
        """
        The find method must find a helix with a repeat count of two when it
        is not at the start of the sequence.
        """
        read = AARead('id', 'AAAFRRFRRF')
        landmark = AlphaHelix_3_10()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('AlphaHelix_3_10', 'B', 3, 7, 2)], result)

    def testHelixTwoRepeatsWithNonZeroOffsetWithSuffix(self):
        """
        The find method must find a helix with a repeat count of two when it
        is not at the start of the sequence and there is also a non-matching
        suffix.
        """
        read = AARead('id', 'AAAFRRFRRFAAA')
        landmark = AlphaHelix_3_10()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('AlphaHelix_3_10', 'B', 3, 7, 2)], result)

    def testHelixThreeRepeats(self):
        """
        The find method must find a helix with a repeat count of three.
        """
        read = AARead('id', 'FRRFRRFRRF')
        landmark = AlphaHelix_3_10()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('AlphaHelix_3_10', 'B', 0, 10, 3)], result)

    def testTwoHelices(self):
        """
        The find method must find more than one helix.
        """
        read = AARead('id', 'FRRFRRFRFRFRFRFRFRFRFRFRFRFRFFRRFRRFRRF')
        landmark = AlphaHelix_3_10()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('AlphaHelix_3_10', 'B', 0, 7, 2),
                          Landmark('AlphaHelix_3_10', 'B', 29, 10, 3)],
                         result)
