from unittest import TestCase

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
        The find method must find one motif.
        """
        read = AARead('id', 'EGGELGYAAAAAAAA')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS60002', 0, 7)], result)

    def testFindOneMotifMiddle(self):
        """
        The find method must find one motif.
        """
        read = AARead('id', 'AAAAAAAAEGGELGYAAAAAAAA')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS60002', 8, 7)], result)

    def testFindOneMotifEnd(self):
        """
        The find method must find one motif.
        """
        read = AARead('id', 'AAAAAAAAEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS60002', 8, 7)], result)

    def testFindTwoIdenticalMotifsNotAdjacent(self):
        """
        The find method must find two identical motifs.
        """
        read = AARead('id', 'EGGELGYAAAAAAAAEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS60002', 0, 7),
                          Landmark('Prosite', 'PS60002', 15, 7)], result)

    def testFindTwoIdenticalMotifsAdjacent(self):
        """
        The find method must find two identical motifs.
        """
        read = AARead('id', 'EGGELGYEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS60002', 0, 7),
                          Landmark('Prosite', 'PS60002', 7, 7)], result)

    def testFindTwoDifferentMotifs(self):
        """
        The find method must find two different motifs.
        """
        read = AARead('id', 'RGDFRRRFEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS00016', 0, 3),
                          Landmark('Prosite', 'PS60002', 8, 7)], result)

    def testTwoDifferentMotifsTwoOccurences(self):
        """
        The find method must find two different motifs with two occurences.
        """
        read = AARead('id', 'EGGELGYARGDAAARGDAAEGGELGY')
        landmark = Prosite()
        result = list(landmark.find(read))
        self.assertEqual([Landmark('Prosite', 'PS00016', 8, 3),
                          Landmark('Prosite', 'PS00016', 14, 3),
                          Landmark('Prosite', 'PS60002', 0, 7),
                          Landmark('Prosite', 'PS60002', 19, 7)], result)
