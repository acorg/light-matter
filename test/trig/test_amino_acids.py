from unittest import TestCase

from dark.reads import AARead

from light.features import TrigPoint
from light.trig.amino_acids import AminoAcids


class TestAminoAcids(TestCase):
    """
    Tests for the TrigPoint.AminoAcids class.
    """
    def testFindWithoutAA(self):
        """
        If the desired aa is not present, find nothing
        """
        read = AARead('id', 'ASAAAAASAAA')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([], result)

    def testFindWithOneAA(self):
        """
        If the desired aa is not present, find nothing
        """
        read = AARead('id', 'ASAACAASAAA')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([TrigPoint('AminoAcids', 'M', 4)], result)

    def testFindWithTwoAA(self):
        """
        If the desired aa is not present, find nothing
        """
        read = AARead('id', 'ASAACAACAAA')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([TrigPoint('AminoAcids', 'M', 4),
                          TrigPoint('AminoAcids', 'M', 7)], result)

    def testFindWithTwoAABeginningAndEnd(self):
        """
        If the desired aa is not present, find nothing
        """
        read = AARead('id', 'CASAAAAAAAC')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([TrigPoint('AminoAcids', 'M', 0),
                          TrigPoint('AminoAcids', 'M', 10)], result)
