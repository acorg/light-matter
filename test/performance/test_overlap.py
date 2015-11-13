import builtins
from unittest import TestCase
from unittest.mock import patch

from light.performance.overlap import SSAARead, CalculateOverlap

from dark.reads import Reads

from ..mocking import mockOpen


class TestSSAARead(TestCase):
    """
    Tests for the SSAARead class.
    """
    def testSSAAReadCorrectAttributes(self):
        """
        An SSAARead must have the correct attributes.
        """
        read = SSAARead('id', 'AFGGCED', 'HHH  HH')
        self.assertEqual('id', read.id)
        self.assertEqual('AFGGCED', read.sequence)
        self.assertEqual('HHH  HH', read.structure)

    def testReadsWithSSAAReads(self):
        """
        It must be possible to make a dark.Reads object out of SSAAReads
        with the right length.
        """
        reads = Reads()
        read1 = SSAARead('id1', 'AFGGCED', 'HHHHHHH')
        read2 = SSAARead('id2', 'AFGGKLL', 'HHHHIII')
        reads.add(read1)
        reads.add(read2)
        self.assertEqual(2, len(reads))


class TestCalculateOverlap(TestCase):
    """
    Tests for the CalculateOverlap class.
    """
    def testEmptyFile(self):
        """
        An empty pdb file must lead to an empty C{dark.reads.Reads} object.
        """
        reads = Reads()
        mockOpener = mockOpen()
        with patch.object(builtins, 'open', mockOpener):
            co = CalculateOverlap('pdb.fasta')
            self.assertEqual(list(reads), list(co.SSAAReads))

    def xtestOneReadFileParsing(self):
        """
        Parsing of the pdb file with one read must work correctly.
        """
        reads = Reads()
        read = SSAARead('seq1', 'REDD', 'HH--')
        reads.add(read)
        data = '\n'.join(['>seq1', 'REDD', '>str1', 'HH--'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            co = CalculateOverlap('pdb.fasta')
            self.assertEqual(1, len(co.SSAAReads))
            self.assertEqual(list(reads)[0].id, list(co.SSAAReads)[0].id)
            self.assertEqual(list(reads)[0].sequence,
                             list(co.SSAAReads)[0].sequence)
            self.assertEqual(list(reads)[0].structure,
                             list(co.SSAAReads)[0].structure)

    def testTwoReadsFileParsing(self):
        """
        Parsing of the pdb file with two reads must work correctly.
        """
        reads = Reads()
        read1 = SSAARead('seq1', 'REDD', 'HH--')
        read2 = SSAARead('seq2', 'REAA', 'HHEE')
        reads.add(read1)
        reads.add(read2)
        data = '\n'.join(['>seq1', 'REDD', '>str1', 'HH--',
                         '>seq2', 'REAA', '>str2', 'HHEE'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            co = CalculateOverlap('pdb.fasta')
            self.assertEqual(2, len(co.SSAAReads))
            self.assertEqual(list(reads)[0].id, list(co.SSAAReads)[0].id)
            self.assertEqual(list(reads)[0].sequence,
                             list(co.SSAAReads)[0].sequence)
            self.assertEqual(list(reads)[0].structure,
                             list(co.SSAAReads)[0].structure)
            self.assertEqual(list(reads)[1].id, list(co.SSAAReads)[1].id)
            self.assertEqual(list(reads)[1].sequence,
                             list(co.SSAAReads)[1].sequence)
            self.assertEqual(list(reads)[1].structure,
                             list(co.SSAAReads)[1].structure)

    def testGetFeatures(self):
        """
        Calling getFeatures on an SSAARead must give the correct result.
        """
        sequence = 'SMEQVAMELRLTELTRLLRSVLDQLQDKDPARIFAQPVSLKEVPDYLDHIKHPMD'
        structure = '-HHHHHHHHHHHHHHHHHHHHHHHHHHT-TT-TTSS---TTT-TTHHHH-SS---'

        ssAARead = SSAARead('5AMF', sequence, structure)
        aa, ss, int_ = CalculateOverlap.getFeatures(ssAARead)
        self.assertEqual({
            'AlphaHelix_pi': set(),
            'Betastrand': set(),
            'AlphaHelix': set(),
            'GOR4AlphaHelix': {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                               18, 19, 20, 21, 22, 23, 24, 29, 30, 31, 32, 33},
            'AlphaHelix_3_10': {21, 22, 23, 24, 25, 26, 27, 28, 29, 30},
            'GOR4BetaStrand': {53, 54}}, aa)
        self.assertEqual({
            '-': {0, 36, 37, 38, 42, 49, 52, 53, 54, 28, 31},
            'H': {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                  18, 19, 20, 21, 22, 23, 24, 25, 26, 45, 46, 47, 48},
            'I': set(),
            'G': set(),
            'E': set(),
            'S': {51, 34, 35, 50},
            'T': {32, 33, 39, 40, 41, 43, 44, 27, 29, 30}}, ss)
        self.assertEqual({
            'GOR4AlphaHelix_HGI': {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                   14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                   25, 26, 29, 30, 31, 32, 33, 45, 46, 47, 48},
            'Betastrand_E': set(),
            'AlphaHelix_H': {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                             16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 45,
                             46, 47, 48},
            'GOR4BetaStrand_E': {53, 54},
            'AlphaHelix_3_10_G': {21, 22, 23, 24, 25, 26, 27, 28, 29, 30},
            'AlphaHelix_pi_I': set()}, int_)

    def testCorrectFractionsMustBeCalculated(self):
        """
        Test that the correct fractions are calculated.
        """
        sequence = 'SMEQVAMELRLTELTRLLRSVLDQLQDKDPARIFAQPVSLKEVPDYLDHIKHPMD'
        structure = ' HHHHHHHHHHHHHHHHHHHHHHHHHHT TT TTSS   TTT TTHHHH SS   '

        ssAARead = SSAARead('5AMF', sequence, structure)
        aa, ss, int_ = CalculateOverlap.getFeatures(ssAARead)
        (alphaHelix, alphaHelix_pi, alphaHelix_3_10, gor4AlphaHelix,
         betaStrand, gor4BetaStrand) = CalculateOverlap.calculateFraction(aa,
                                                                          ss,
                                                                          int_)
        self.assertEqual(0.0, alphaHelix)
        self.assertEqual(0.0, alphaHelix_pi)
        self.assertEqual(1.0, alphaHelix_3_10)
        self.assertAlmostEqual(0.7142857, gor4AlphaHelix)
        self.assertEqual(0.0, betaStrand)
        self.assertEqual(1.0, gor4BetaStrand)
