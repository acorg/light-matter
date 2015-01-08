from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.turn import Turn


class TestTurn(TestCase):
    """
    Tests for the Landmark.Turn class.
    """
    def testClassAttributes(self):
        """
        The Turn class attributes must be as expected.
        """
        self.assertEqual('Turn', Turn.NAME)
        self.assertEqual('T', Turn.SYMBOL)

    def testNoTurn(self):
        """
        The find method must return an empty generator when no turn is present.
        """
        read = AARead('id', 'PAPAPA')
        landmark = Turn()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testAtStartOfSequence(self):
        """
        The find method must find a turn that begins at the start of the
        sequence.
        """
        read = AARead('id', 'NPGYAAAAAA')
        landmark = Turn()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Turn', 'T', 0, 4)], result)

    def testAtEndOfSequence(self):
        """
        The find method must find turn that occurs at the end of the sequence.
        """
        read = AARead('id', 'AAAADPDG')
        landmark = Turn()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Turn', 'T', 4, 4)], result)

    def testInMiddleOfSequence(self):
        """
        The find method must find a turn that starts in the middle of the
        sequence.
        """
        read = AARead('id', 'AADPDGAA')
        landmark = Turn()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Turn', 'T', 2, 4)], result)

    def testMultipleMatches(self):
        """
        The find method must find several turn sequences in the same read.
        """
        read = AARead('id', 'NPNWAACSDYAADKAY')
        landmark = Turn()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Turn', 'T', 0, 4),
                          Landmark('Turn', 'T', 6, 4)],
                         result)
