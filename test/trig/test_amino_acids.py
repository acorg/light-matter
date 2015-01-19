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
        If the desired amino acid is not present, return an empty list.
        """
        read = AARead('id', 'ASAAAAASAAA')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([], result)

    def testFindWithOneAAW(self):
        """
        If one of the desired amino acids (W) is present, return one trig
        point.
        """
        read = AARead('id', 'ASAAWAASAAA')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([TrigPoint('AminoAcids', 'M', 4)], result)

    def testFindWithTwoAAW(self):
        """
        If two of the desired amino acids (W) are present, return two trig
        points.
        """
        read = AARead('id', 'ASAAWAAWAAA')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([TrigPoint('AminoAcids', 'M', 4),
                          TrigPoint('AminoAcids', 'M', 7)], result)

    def testFindWithTwoAABeginningAndEndW(self):
        """
        If two of the desired amino acids are present (W) at the beginning and
        end of the sequence, return two trig points.
        """
        read = AARead('id', 'WASAAAAAAAW')
        aas = AminoAcids()
        result = list(aas.find(read))
        self.assertEqual([TrigPoint('AminoAcids', 'M', 0),
                          TrigPoint('AminoAcids', 'M', 10)], result)

    def testFindRightNonDefaultAA(self):
        """
        If the default amino acid list is changed, return the right trig point.
        """
        read = AARead('id', 'ASAACAAWAAA')
        aas = AminoAcids()
        result = list(aas.find(read, aa=['S']))
        self.assertEqual([TrigPoint('AminoAcids', 'M', 1)], result)
