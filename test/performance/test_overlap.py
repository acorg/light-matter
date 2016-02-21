import six
from six.moves import builtins
from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

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

    def testOneReadFileParsing(self):
        """
        Parsing of the pdb file with one read must work correctly.
        """
        read = SSAARead('seq1', 'REDD', 'HH--')
        data = '\n'.join(['>seq1', 'REDD', '>str1', 'HH--'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            co = CalculateOverlap('pdb.fasta')
            self.assertEqual(1, len(co.SSAAReads))
            self.assertEqual(read.id, list(co.SSAAReads)[0].id)
            self.assertEqual(read.sequence,
                             list(co.SSAAReads)[0].sequence)
            self.assertEqual(read.structure,
                             list(co.SSAAReads)[0].structure)

    def testOddNumberOfRecords(self):
        """
        Trying to parse a PDB file with an odd number of records must raise a
        ValueError.
        """
        data = '\n'.join(['>seq1', 'REDD',
                          '>str1', 'HH--',
                          '>seq2', 'REAA'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = "^Structure file 'x.fasta' has an odd number of records\.$"
            six.assertRaisesRegex(self, ValueError, error, CalculateOverlap,
                                  'x.fasta')

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
            self.assertEqual(read1.id, list(co.SSAAReads)[0].id)
            self.assertEqual(read1.sequence,
                             list(co.SSAAReads)[0].sequence)
            self.assertEqual(read1.structure,
                             list(co.SSAAReads)[0].structure)
            self.assertEqual(read2.id, list(co.SSAAReads)[1].id)
            self.assertEqual(read2.sequence,
                             list(co.SSAAReads)[1].sequence)
            self.assertEqual(read2.structure,
                             list(co.SSAAReads)[1].structure)

    def testGetFeatures(self):
        """
        Calling getFeatures on an SSAARead must give the correct result.
        """
        sequence = 'SMEQVAMELRLTELTRLLRSVLDQLQDKDPARIFAQPVSLKEVPDYLDHIKHPMD'
        structure = '-HHHHHHHHHHHHHHHHHHHHHHHHHHT-TT-TTSS---TTT-TTHHHH-SS---'

        ssAARead = SSAARead('5AMF', sequence, structure)
        seqF, commons, totals = CalculateOverlap.getFeatures(ssAARead)
        self.assertEqual({
            'AlphaHelix': set(),
            'AlphaHelix_3_10': {21, 22, 23, 24, 25, 26, 27, 28, 29, 30},
            'AlphaHelix_pi': set(),
            'AminoAcidsLm': set(),
            'BetaStrand': set(),
            'BetaTurn': set(),
            'GOR4AlphaHelix': {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                               18, 19, 20, 21, 22, 23, 24, 29, 30, 31, 32, 33},
            'GOR4BetaStrand': {53, 54},
            'GOR4Coil': {0, 1, 2, 3, 4, 25, 26, 27, 28, 34, 35, 36, 37, 38, 39,
                         40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52},
            'Prosite': {32, 38, 39, 40, 41, 14, 15, 16, 19, 20, 21, 22, 30,
                        31},
            'AminoAcids': set(),
            'IndividualPeaks': {9, 31},
            'IndividualTroughs': {49, 10},
            'Peaks': {3, 35, 5, 38, 7, 40, 9, 12, 44, 15, 47, 18, 50, 22, 27,
                      31},
            'Troughs': {1, 33, 4, 37, 6, 39, 8, 10, 42, 13, 46, 49, 21, 53, 24,
                        30},
            '-': {0, 36, 37, 38, 42, 49, 52, 53, 54, 28, 31},
            'E': set(),
            'G': set(),
            'H': {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                  18, 19, 20, 21, 22, 23, 24, 25, 26, 45, 46, 47, 48},
            'I': set(),
            'S': {51, 34, 35, 50},
            'T': {32, 33, 39, 40, 41, 43, 44, 27, 29, 30}}, seqF)
        self.assertEqual({49, 10}, commons['IndividualTroughs-Troughs'])
        self.assertEqual({1, 33, 4, 37, 6, 39, 8, 10, 42, 13, 46, 49, 21, 53,
                          24, 30}, totals['IndividualTroughs-Troughs'])
        self.assertEqual(190, len(commons))
        self.assertEqual(190, len(totals))
