import builtins
from unittest import TestCase
from unittest.mock import patch
from test.mocking import mockOpen

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.prosite import Prosite


class TestProsite(TestCase):
    """
    Tests for the Landmark.Prosite class.
    """

    # A simple prosite JSON db with three patterns, to save on reading the
    # default large prosite file on every test.
    SIMPLE_DB = (
        '{"accession": "00016", "pattern": "RGD"}'
        '\n'
        '{"accession": "60002", "pattern": "EGGELGY"}'
        '\n'
        '{"accession": "00021", "pattern": "[FY]{1,2}C[RH].{3,4}[^WY]C"}'
    )

    def testFindWithoutMotif(self):
        """
        The find method must return an empty generator when no motif is
        present.
        """
        read = AARead('id', 'FFFFFFFFFFFFFFFFFFFFFFFFFFF')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch.object(builtins, 'open', mockOpener):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([], result)

    def testFindOneMotifBeginning(self):
        """
        The find method must find one motif at the beginning of a sequence.
        """
        read = AARead('id', 'EGGELGYAAAAAAAA')
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch.object(builtins, 'open', mockOpener):
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
        with patch.object(builtins, 'open', mockOpener):
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
        with patch.object(builtins, 'open', mockOpener):
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
        with patch.object(builtins, 'open', mockOpener):
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
        with patch.object(builtins, 'open', mockOpener):
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
        with patch.object(builtins, 'open', mockOpener):
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
        with patch.object(builtins, 'open', mockOpener):
            landmark = Prosite()
            result = list(landmark.find(read))
            self.assertEqual([Landmark('Prosite', 'PS', 8, 3, '00016'),
                              Landmark('Prosite', 'PS', 14, 3, '00016'),
                              Landmark('Prosite', 'PS', 0, 7, '60002'),
                              Landmark('Prosite', 'PS', 19, 7, '60002')],
                             result)

    def testFindWithPattern(self):
        """
        The find method must work on a regex that contains alternation (e.g.,
        [ABC]), negated alternation (e.g., [^ABC]), and repetition (e.g.,
        .(7,8) or [ABC](1,3)).
        """
        # The pattern we're matching is [FY](1,2)C[RH].(3,4)[^WY]C
        mockOpener = mockOpen(read_data=self.SIMPLE_DB)
        with patch.object(builtins, 'open', mockOpener):
            landmark = Prosite()
            # Put 2 X's in for .(3,4)
            result = list(landmark.find(AARead('id', 'FCRXXBC')))
            self.assertEqual([], result)
            # Put 3 X's in for .(3,4)
            result = list(landmark.find(AARead('id', 'FCRXXXBC')))
            self.assertEqual([Landmark('Prosite', 'PS', 0, 8, '00021')],
                             result)
            # Put 4 X's in for .(3,4)
            result = list(landmark.find(AARead('id', 'FCRXXXXBC')))
            self.assertEqual([Landmark('Prosite', 'PS', 0, 9, '00021')],
                             result)
            # Put 5 X's in for .(3,4)
            result = list(landmark.find(AARead('id', 'FCRXXXXXBC')))
            self.assertEqual([], result)
            # Put X in for [FY]
            result = list(landmark.find(AARead('id', 'XCRXXXBC')))
            self.assertEqual([], result)
            # Put W in for [^WY]
            result = list(landmark.find(AARead('id', 'FCRXXXXWC')))
            self.assertEqual([], result)
            # Put YY in for [FY](1,2)
            result = list(landmark.find(AARead('id', 'YYCRXXXBC')))
            self.assertEqual([Landmark('Prosite', 'PS', 0, 9, '00021')],
                             result)

    def testDefaultDatabase(self):
        """
        When no database file is passed to Prosite, it must find and open
        the default database, which must have the correct length.
        """
        landmark = Prosite()
        # The 1309 is based on Prosite database 20.119 of 14-Oct-2015
        self.assertEqual(1309, len(landmark.database))
