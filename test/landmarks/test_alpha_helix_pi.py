from unittest import TestCase

from dark.reads import AARead

from light.landmark import Landmark
from light.landmarks.alpha_helix_pi import AlphaHelix_pi


class TestAlphaHelix_pi(TestCase):
    """
    Tests for the Landmark.AlphaHelix_pi class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right properties
        string.
        """
        landmark = AlphaHelix_pi()
        result = landmark.convertAAToHydrophobicHydrophilic('ASDGEAHSDTDSCV')
        self.assertEqual('OIIIIOOIIOIIOO', result)

    def testFindWithoutHelix(self):
        """
        The find method must return an empty generator when no helix is
        present.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        landmark = AlphaHelix_pi()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testFindWithOneRepeat(self):
        """
        The find method must return an empty generator if there is only one
        instance of the OII pattern
        """
        read = AARead('id', 'FRRRR')
        landmark = AlphaHelix_pi()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testHelixTwoRepeats(self):
        """
        The find method must find a helix with a repeat count of two.
        """
        read = AARead('id', 'FRRRRFRRRRF')
        landmark = AlphaHelix_pi()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('C', 0, 11, 2)], result)

    def testHelixTwoRepeatsWithNonZeroOffset(self):
        """
        The find method must find a helix with a repeat count of two when it
        is not at the start of the sequence.
        """
        read = AARead('id', 'AAAFRRRRFRRRRF')
        landmark = AlphaHelix_pi()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('C', 3, 11, 2)], result)

    def testHelixTwoRepeatsWithNonZeroOffsetWithSuffix(self):
        """
        The find method must find a helix with a repeat count of two when it
        is not at the start of the sequence and there is also a non-matching
        suffix.
        """
        read = AARead('id', 'AAAFRRRRFRRRRFAAA')
        landmark = AlphaHelix_pi()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('C', 3, 11, 2)], result)

    def testHelixThreeRepeats(self):
        """
        The find method must find a helix with a repeat count of three.
        """
        read = AARead('id', 'FRRRRFRRRRFRRRRF')
        landmark = AlphaHelix_pi()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('C', 0, 16, 3)], result)

    def testTwoHelices(self):
        """
        The find method must find more than one helix.
        """
        read = AARead('id', 'FRRRRFRRRRFRFRFRFRFRFRFRFRFRFRFFRRRRFRRRRFRRRRF')
        landmark = AlphaHelix_pi()
        result = list(landmark.find(read))
        self.assertEqual(
            [Landmark('C', 0, 11, 2), Landmark('C', 31, 16, 3)], result)
