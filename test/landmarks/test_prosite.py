from unittest import TestCase
from test.mocking import mockOpen
from mock import patch

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.prosite import Prosite


class TestProsite(TestCase):
    """
    Tests for the Landmark.Prosite class.
    """

    # A simple prosite JSON db with two patterns, to save on reading the
    # default large prosite file on every test.
    SIMPLE_DB = ('{"accession": "PS00016", "pattern": "RGD"}'
                 '\n'
                 '{"accession": "PS60002", "pattern": "EGGELGY"}')

    def testFindWithoutMotif(self):
        """
        The find method must return an empty generator when no motif is
        present.
        """
        read = AARead('id', 'FFFFFFFFFFFFFFFFFFFFFFFFFFF')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch('__builtin__.open', mockOpener, create=True):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([], result)

    def testFindOneMotifBeginning(self):
        """
        The find method must find one motif at the beginning of a sequence.
        """
        read = AARead('id', 'EGGELGYAAAAAAAA')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch('__builtin__.open', mockOpener, create=True):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([Landmark('Prosite', 'PS', 0, 7, '60002')],
                             result)

    def testFindOneMotifMiddle(self):
        """
        The find method must find one motif in the middle of a sequence.
        """
        read = AARead('id', 'AAAAAAAAEGGELGYAAAAAAAA')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch('__builtin__.open', mockOpener, create=True):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([Landmark('Prosite', 'PS', 8, 7, '60002')],
                             result)

    def testFindOneMotifEnd(self):
        """
        The find method must find one motif at the end of a sequence.
        """
        read = AARead('id', 'AAAAAAAAEGGELGY')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch('__builtin__.open', mockOpener, create=True):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([Landmark('Prosite', 'PS', 8, 7, '60002')],
                             result)

    def testFindTwoIdenticalMotifsNotAdjacent(self):
        """
        The find method must find two identical motifs.
        """
        read = AARead('id', 'EGGELGYAAAAAAAAEGGELGY')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch('__builtin__.open', mockOpener, create=True):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([Landmark('Prosite', 'PS', 0, 7, '60002'),
                              Landmark('Prosite', 'PS', 15, 7, '60002')],
                             result)

    def testFindTwoIdenticalMotifsAdjacent(self):
        """
        The find method must find two identical motifs.
        """
        read = AARead('id', 'EGGELGYEGGELGY')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch('__builtin__.open', mockOpener, create=True):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([Landmark('Prosite', 'PS', 0, 7, '60002'),
                              Landmark('Prosite', 'PS', 7, 7, '60002')],
                             result)

    def testFindTwoDifferentMotifs(self):
        """
        The find method must find two different motifs.
        """
        read = AARead('id', 'RGDFRRRFEGGELGY')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch('__builtin__.open', mockOpener, create=True):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([Landmark('Prosite', 'PS', 0, 3, '00016'),
                              Landmark('Prosite', 'PS', 8, 7, '60002')],
                             result)

    def testTwoDifferentMotifsTwoOccurences(self):
        """
        The find method must find two different motifs with two occurences.
        """
        read = AARead('id', 'EGGELGYARGDAAARGDAAEGGELGY')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch('__builtin__.open', mockOpener, create=True):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([Landmark('Prosite', 'PS', 8, 3, '00016'),
                              Landmark('Prosite', 'PS', 14, 3, '00016'),
                              Landmark('Prosite', 'PS', 0, 7, '60002'),
                              Landmark('Prosite', 'PS', 19, 7, '60002')],
                             result)

    def testAccessionPrefix(self):
        """
        The find method must raise an error if the first two letters of the
        accession are not 'PS'.
        """
        data = '{"accession": "XX00016", "pattern": "EGGELGY"}'
        mockOpener = mockOpen(read_data=data)
        with patch('__builtin__.open', mockOpener, create=True):
            self.assertRaises(AssertionError, Prosite, 'database.dat')

    def testDuplicateAccession(self):
        """
        The find method must raise AssertionError if the database contains a
        duplicate accession string.
        """
        data = ('{"accession": "XX00016", "pattern": "EGGELGY"}'
                '\n'
                '{"accession": "XX00016", "pattern": "EGGELGY"}')
        mockOpener = mockOpen(read_data=data)
        with patch('__builtin__.open', mockOpener, create=True):
            self.assertRaises(AssertionError, Prosite, 'database.dat')

    def testDefaultDatabase(self):
        """
        When no database file is passed to Prosite, it must find and open
        the default database, which must have the correct length.
        """
        landmark = Prosite()
        self.assertEqual(1308, len(landmark.database))
