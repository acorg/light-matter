from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.beta_turn import BetaTurn


class TestBetaTurn(TestCase):
    """
    Tests for the Landmark.BetaTurn class.
    """
    def testClassAttributes(self):
        """
        The BetaTurn class attributes must be as expected.
        """
        self.assertEqual('BetaTurn', BetaTurn.NAME)
        self.assertEqual('BT', BetaTurn.SYMBOL)

    def testNoBetaTurn(self):
        """
        The find method must return an empty generator when no beta turn is
        present.
        """
        read = AARead('id', 'PAPAPA')
        landmark = BetaTurn()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testAtStartOfSequence(self):
        """
        The find method must find a beta turn that begins at the start of the
        sequence.
        """
        read = AARead('id', 'NPGYAAAAAA')
        landmark = BetaTurn()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('BetaTurn', 'BT', 0, 4)], result)

    def testAtEndOfSequence(self):
        """
        The find method must find a beta turn that occurs at the end of the
        sequence.
        """
        read = AARead('id', 'AAAADPDG')
        landmark = BetaTurn()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('BetaTurn', 'BT', 4, 4)], result)

    def testInMiddleOfSequence(self):
        """
        The find method must find a beta turn that starts in the middle of the
        sequence.
        """
        read = AARead('id', 'AADPDGAA')
        landmark = BetaTurn()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('BetaTurn', 'BT', 2, 4)], result)

    def testMultipleMatches(self):
        """
        The find method must find several beta turn sequences in the same read.
        """
        read = AARead('id', 'NPNWAACSDYAADKAY')
        landmark = BetaTurn()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('BetaTurn', 'BT', 0, 4),
                          Landmark('BetaTurn', 'BT', 6, 4)],
                         result)

    def testNoLongTurn(self):
        """
        The returned beta turn must not be longer than four residues.
        """
        read = AARead('id', 'NPGYNPAAAAAA')
        landmark = BetaTurn()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('BetaTurn', 'BT', 0, 4)], result)
