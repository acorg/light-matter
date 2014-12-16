from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.beta_strand import BetaStrand


class TestBetaStrand(TestCase):
    """
    Tests for the Landmark.BetaStrand class.
    """
    def testClassAttributes(self):
        """
        The BetaStrand class attributes must be as expected.
        """
        self.assertEqual('BetaStrand', BetaStrand.NAME)
        self.assertEqual('S', BetaStrand.SYMBOL)
        self.assertEqual('VILMCFYWQTHKR', BetaStrand.BETA_STRAND_AAs)
        self.assertEqual(3, BetaStrand.MIN_LENGTH)

    def testNoBetaStrand(self):
        """
        The find method must return an empty generator when no beta strand is
        present.
        """
        read = AARead('id', 'PAPAPA')
        landmark = BetaStrand()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testBelowMinimumLength(self):
        """
        The find method must not find a too-short (<3 AA) sequence of beta
        strand amino acids.
        """
        read = AARead('id', 'VI')
        landmark = BetaStrand()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testMinimumLength(self):
        """
        The find method must find a minimal length (3) sequence of beta strand
        amino acids.
        """
        read = AARead('id', 'VIL')
        landmark = BetaStrand()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('BetaStrand', 'S', 0, 3, 3)], result)

    def testNotAtStartOfSequence(self):
        """
        The find method must find a sequence of beta strand amino acids that
        does not begin at the start of the sequence.
        """
        read = AARead('id', 'PVILVIL')
        landmark = BetaStrand()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('BetaStrand', 'S', 1, 6, 6)], result)

    def testInMiddleOfSequence(self):
        """
        The find method must find a sequence of beta strand amino acids that
        starts in the middle of the sequence.
        """
        read = AARead('id', 'PAVILMCFYPA')
        landmark = BetaStrand()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('BetaStrand', 'S', 2, 7, 7)], result)

    def testMultipleMatches(self):
        """
        The find method must find several beta strand amino acids sequences
        in the same read.
        """
        read = AARead('id', 'PAVILPAPAVILVILPA')
        landmark = BetaStrand()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('BetaStrand', 'S', 2, 3, 3),
                          Landmark('BetaStrand', 'S', 9, 6, 6)],
                         result)
