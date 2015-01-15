from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.amino_acids import AminoAcids


class TestAminoAcids(TestCase):
    """
    Tests for the TrigPoint.AminoAcids class.
    """
    def testFindWithoutAA(self):
        """
        If the desired amino acid is not present, return an empty list.
        """
        read = AARead('id', 'ASAAAAASAAA')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([], result)

    def testFindWithOneAAC(self):
        """
        If one of the desired amino acids (C) is present, return one landmark.
        """
        read = AARead('id', 'ASAACAASAAA')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([Landmark('AminoAcids', 'N', 4, 1)], result)

    def testFindWithTwoAAC(self):
        """
        If two of the desired amino acids (C) are present, return two
        landmarks.
        """
        read = AARead('id', 'ASAACAACAAA')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([Landmark('AminoAcids', 'N', 4, 1),
                          Landmark('AminoAcids', 'N', 7, 1)], result)

    def testFindWithTwoAABeginningAndEndC(self):
        """
        If two of the desired amino acids are present (C) at the beginning and
        end of the sequence, return two landmarks.
        """
        read = AARead('id', 'CASAAAAAAAC')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([Landmark('AminoAcids', 'N', 0, 1),
                          Landmark('AminoAcids', 'N', 10, 1)], result)

    def testFindRightNonDefaultAA(self):
        """
        If the default amino acid list is changed, return the right trig point.
        """
        read = AARead('id', 'ASAACAAWAAA')
        aas = AminoAcids()
        result = list(aas.find(read, aa=['S']))
        self.assertEqual([Landmark('AminoAcids', 'N', 1, 1)], result)
