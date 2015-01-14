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
    def testFindWithoutMotif(self):
        """
        The find method must return an empty generator when no motif is
        present.
        """
        read = AARead('id', 'FFFFFFFFFFFFFFFFFFFFFFFFFFF')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testFindOneMotifBeginning(self):
        """
        The find method must find one motif at the beginning of a sequence.
        """
        read = AARead('id', 'EGGELGYAAAAAAAA')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 7, '60002')], result)

    def testFindOneMotifMiddle(self):
        """
        The find method must find one motif in the middle of a sequence.
        """
        read = AARead('id', 'AAAAAAAAEGGELGYAAAAAAAA')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 8, 7, '60002')], result)

    def testFindOneMotifEnd(self):
        """
        The find method must find one motif at the end of a sequence.
        """
        read = AARead('id', 'AAAAAAAAEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 8, 7, '60002')], result)

    def testFindTwoIdenticalMotifsNotAdjacent(self):
        """
        The find method must find two identical motifs.
        """
        read = AARead('id', 'EGGELGYAAAAAAAAEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 7, '60002'),
                          Landmark('Prosite', 'PS', 15, 7, '60002')], result)

    def testFindTwoIdenticalMotifsAdjacent(self):
        """
        The find method must find two identical motifs.
        """
        read = AARead('id', 'EGGELGYEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 7, '60002'),
                          Landmark('Prosite', 'PS', 7, 7, '60002')], result)

    def testFindTwoDifferentMotifs(self):
        """
        The find method must find two different motifs.
        """
        read = AARead('id', 'RGDFRRRFEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 3, '00016'),
                          Landmark('Prosite', 'PS', 8, 7, '60002')], result)

    def testTwoDifferentMotifsTwoOccurences(self):
        """
        The find method must find two different motifs with two occurences.
        """
        read = AARead('id', 'EGGELGYARGDAAARGDAAEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 8, 3, '00016'),
                          Landmark('Prosite', 'PS', 14, 3, '00016'),
                          Landmark('Prosite', 'PS', 0, 7, '60002'),
                          Landmark('Prosite', 'PS', 19, 7, '60002')], result)

    def testRaise(self):
        """
        The find method must raise an error if the first two letters of the
        accession are not 'PS'.
        """
        read = AARead('id', 'EGGELGYA')
        data = '{"accession": "XX00016", "pattern": "EGGELGY"}'
        mockOpener = mockOpen(read_data=data)
        with patch('__builtin__.open', mockOpener, create=True):
            landmark = Prosite(databaseFile='database.dat')
            self.assertRaises(AssertionError, landmark.find, read)
