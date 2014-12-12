from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.individual_troughs import IndividualTroughs


class TestIndividualTroughs(TestCase):
    """
    Tests for the TrigPoint.IndividualTroughs class.
    """
    def testFindWithoutIndividualTrough(self):
        """
        The find method must return an empty generator when no individual
        trough is present.
        """
        read = AARead('id', 'RRRRRRRRRRRRRRR')
        troughs = IndividualTroughs()
        result = list(troughs.find(read))
        self.assertEqual([], result)

    def testTwoIndividualTroughs(self):
        """
        The find method must find two individual troughs.
        """
        read = AARead('id', 'TMTTTTTMTTT')
        troughs = IndividualTroughs()
        result = list(troughs.find(read))
        self.assertEqual([TrigPoint('IndividualTroughs', 'J', 1),
                          TrigPoint('IndividualTroughs', 'J', 7)], result)

    def testNoIndividualTroughAtBeginning(self):
        """
        The find method must not find an individual trough at the beginning of
        the sequence.
        """
        read = AARead('id', 'VAAAAAAAAA')
        troughs = IndividualTroughs()
        result = list(troughs.find(read))
        self.assertEqual([], result)

    def testNoIndividualTroughAtEnd(self):
        """
        The find method must not find an individual trough at the end of the
        sequence.
        """
        read = AARead('id', 'AAAAAAV')
        troughs = IndividualTroughs()
        result = list(troughs.find(read))
        self.assertEqual([], result)
