from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.prosite import Prosite, _DATABASE, _loadDatabase


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
        finder = Prosite()
        result = list(finder.find(read))
        self.assertEqual([], result)

    def testFindOneMotifBeginning(self):
        """
        The find method must find one motif at the beginning of a sequence.
        """
        read = AARead('id', 'EGGELGYAAAAAAAA')
        finder = Prosite()
        result = list(finder.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 1, 6, '00008'),
                          Landmark('Prosite', 'PS', 0, 7, '60002')],
                         result)

    def testFindOneMotifMiddle(self):
        """
        The find method must find a motif in the middle of a sequence.
        """
        read = AARead('id', 'AAAAAAAAEGGELGYAAAAAAAA')
        finder = Prosite()
        result = list(finder.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 9, 6, '00008'),
                          Landmark('Prosite', 'PS', 8, 7, '60002')],
                         result)

    def testFindOneMotifEnd(self):
        """
        The find method must find a motif at the end of a sequence.
        """
        read = AARead('id', 'AAAAAAAAEGGELGY')
        finder = Prosite()
        result = list(finder.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 9, 6, '00008'),
                          Landmark('Prosite', 'PS', 8, 7, '60002')],
                         result)

    def testFindTwoIdenticalMotifsNotAdjacent(self):
        """
        The find method must find repeating identical motifs.
        """
        read = AARead('id', 'EGGELGYAAAAAAAAEGGELGY')
        finder = Prosite()
        result = list(finder.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 1, 6, '00008'),
                          Landmark('Prosite', 'PS', 16, 6, '00008'),
                          Landmark('Prosite', 'PS', 0, 7, '60002'),
                          Landmark('Prosite', 'PS', 15, 7, '60002')],
                         result)

    def testFindTwoIdenticalMotifsAdjacent(self):
        """
        The find method must find two identical motifs.
        """
        read = AARead('id', 'EGGELGYEGGELGY')
        finder = Prosite()
        result = list(finder.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 1, 6, '00008'),
                          Landmark('Prosite', 'PS', 8, 6, '00008'),
                          Landmark('Prosite', 'PS', 0, 7, '60002'),
                          Landmark('Prosite', 'PS', 7, 7, '60002')],
                         result)

    def testThreeDifferentMotifsTwoOccurences(self):
        """
        The find method must find three different motifs with two occurences
        each.
        """
        read = AARead('id', 'EGGELGYARGDAAARGDAAEGGELGY')
        finder = Prosite()
        result = list(finder.find(read))
        self.assertEqual([Landmark('Prosite', 'PS', 1, 6, '00008'),
                          Landmark('Prosite', 'PS', 20, 6, '00008'),
                          Landmark('Prosite', 'PS', 8, 3, '00016'),
                          Landmark('Prosite', 'PS', 14, 3, '00016'),
                          Landmark('Prosite', 'PS', 0, 7, '60002'),
                          Landmark('Prosite', 'PS', 19, 7, '60002')],
                         result)

    def testMatchAlteration(self):
        """
        The find method must work on a regex that contains alternation (e.g.,
        [ABC]).
        """
        # Make sure we know what pattern is actually being matched.
        self.assertEqual('[ST].[RK]', _DATABASE[2]['regex'].pattern)
        self.assertEqual('00005', _DATABASE[2]['accession'])

        finder = Prosite()

        result = list(finder.find(AARead('id', 'SXR')))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 3, '00005')], result)

        result = list(finder.find(AARead('id', 'TXR')))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 3, '00005')], result)

        result = list(finder.find(AARead('id', 'SXK')))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 3, '00005')], result)

    def testMatchNegatedAlteration(self):
        """
        The find method must work on a regex that contains negated alternation
        (e.g., [^ABC]).
        """
        # Make sure we know what pattern is actually being matched.
        self.assertEqual('G[^EDRKHPFYW].{2}[STAGCN][^P]',
                         _DATABASE[5]['regex'].pattern)
        self.assertEqual('00008', _DATABASE[5]['accession'])

        finder = Prosite()

        result = list(finder.find(AARead('id', 'GAXXSA')))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 6, '00008')], result)

        result = list(finder.find(AARead('id', 'GEXXSA')))
        self.assertEqual([], result)

        result = list(finder.find(AARead('id', 'GAXXSP')))
        self.assertEqual([], result)

    def testMatchSpecificRepetition(self):
        """
        The find method must work on a regex that contains a pattern with a
        fixed number of repeats (e.g., [RK]{2}.[ST]).
        """
        # Make sure we know what pattern is actually being matched.
        self.assertEqual('[RK]{2}.[ST]',
                         _DATABASE[1]['regex'].pattern)
        self.assertEqual('00004', _DATABASE[1]['accession'])

        finder = Prosite()

        result = list(finder.find(AARead('id', 'RRXS')))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 4, '00004')], result)

        result = list(finder.find(AARead('id', 'RRRRX')))
        self.assertEqual([], result)

    def testMatchRepetitionRange(self):
        """
        The find method must work on a regex that contains a pattern with a
        variable number of repeats (e.g., [RK]{2,4}.[ST]).
        """
        # Make sure we know what pattern is actually being matched.
        self.assertEqual('[RK].{2,3}[DE].{2,3}Y',
                         _DATABASE[4]['regex'].pattern)
        self.assertEqual('00007', _DATABASE[4]['accession'])

        finder = Prosite()

        result = list(finder.find(AARead('id', 'RXXDXXY')))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 7, '00007')], result)

        result = list(finder.find(AARead('id', 'RXDXXY')))
        self.assertEqual([], result)

        result = list(finder.find(AARead('id', 'RXXXDXXY')))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 8, '00007')], result)

        result = list(finder.find(AARead('id', 'RXXXDXXXY')))
        self.assertEqual([Landmark('Prosite', 'PS', 0, 9, '00007')], result)

        result = list(finder.find(AARead('id', 'RXXXDXXXXY')))
        self.assertEqual([], result)

    def testLoadDatabase(self):
        """
        The Prosite database loading function must work and produce the
        expected result.
        """
        db = _loadDatabase()
        # The 1309 is based on Prosite database 20.119 of 14-Oct-2015
        self.assertEqual(1309, len(db))
        expectedKeys = {'regex', 'accession'}
        for motif in db:
            self.assertEqual(expectedKeys, set(motif.keys()))
