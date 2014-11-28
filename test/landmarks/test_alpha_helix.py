from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.alpha_helix import AlphaHelix


class TestAlphaHelix(TestCase):
    """
    Tests for the Landmark.AlphaHelix class.
    """
    def testConversion(self):
        """
        The amino acid sequence must be converted to the right properties
        string.
        """
        landmark = AlphaHelix()
        result = landmark.convertAAToHydrophobicHydrophilic('ASDGEAHSDTDSCV')
        self.assertEqual('OIIIIOOIIOIIOO', result)

    def testFindWithoutHelix(self):
        """
        The find method must return an empty generator when no helix is
        present.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        landmark = AlphaHelix()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testFindWithOneRepeat(self):
        """
        The find method must return an empty generator if there is only one
        instance of the OIII pattern
        """
        read = AARead('id', 'FRRR')
        landmark = AlphaHelix()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testHelixTwoRepeats(self):
        """
        The find method must find a helix with a repeat count of two.
        """
        read = AARead('id', 'FRRRFRRRF')
        landmark = AlphaHelix()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('A', 0, 9, 2)], result)

    def testHelixTwoRepeatsWithNonZeroOffset(self):
        """
        The find method must find a helix with a repeat count of two when it
        is not at the start of the sequence.
        """
        read = AARead('id', 'AAAFRRRFRRRF')
        landmark = AlphaHelix()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('A', 3, 9, 2)], result)

    def testHelixTwoRepeatsWithNonZeroOffsetWithSuffix(self):
        """
        The find method must find a helix with a repeat count of two when it
        is not at the start of the sequence and there is also a non-matching
        suffix.
        """
        read = AARead('id', 'AAAFRRRFRRRFAAA')
        landmark = AlphaHelix()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('A', 3, 9, 2)], result)

    def testHelixThreeRepeats(self):
        """
        The find method must find a helix with a repeat count of three.
        """
        read = AARead('id', 'FRRRFRRRFRRRF')
        landmark = AlphaHelix()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('A', 0, 13, 3)], result)

    def testTwoHelices(self):
        """
        The find method must find more than one helix.
        """
        read = AARead('id', 'FRRRFRRRFRFRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF')
        landmark = AlphaHelix()
        result = list(landmark.find(read))
        self.assertEqual(
            [Landmark('A', 0, 9, 2), Landmark('A', 31, 13, 3)], result)
