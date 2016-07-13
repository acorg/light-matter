from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.elm import EukaryoticLinearMotif, _DATABASE, _loadDatabase


class TestEukaryoticLinearMotif(TestCase):
    """
    Tests for the Landmark.EukaryoticLinearMotif class.
    """

    def testFindWithoutMotif(self):
        """
        The find method must return an empty generator when no motif is
        present.
        """
        read = AARead('id', 'FFFFFFFFFFFFFFFFFFFFFFFFFFF')
        finder = EukaryoticLinearMotif()
        result = list(finder.find(read))
        self.assertEqual([], result)

    def testFindOneMotifBeginning(self):
        """
        The find method must find one motif at the beginning of a sequence.
        """
        read = AARead('id', 'EKENLGYAAAAAAAA')
        finder = EukaryoticLinearMotif()
        result = list(finder.find(read))
        self.assertEqual([Landmark('EukaryoticLinearMotif', 'ELM', 0, 5,
                                   'DEG_APCC_KENBOX_2')], result)

    def testFindOneMotifMiddle(self):
        """
        The find method must find a motif in the middle of a sequence.
        """
        read = AARead('id', 'AAAAAAAAEKENLGYAAAAAAAA')
        finder = EukaryoticLinearMotif()
        result = list(finder.find(read))
        self.assertEqual([Landmark('EukaryoticLinearMotif', 'ELM', 8, 5,
                                   'DEG_APCC_KENBOX_2')], result)

    def testFindOneMotifEnd(self):
        """
        The find method must find a motif at the end of a sequence.
        """
        read = AARead('id', 'AAAAAAAAEGGKENY')
        finder = EukaryoticLinearMotif()
        result = list(finder.find(read))
        self.assertEqual([Landmark('EukaryoticLinearMotif', 'ELM', 10, 5,
                                   'DEG_APCC_KENBOX_2')], result)

    def testFindTwoIdenticalMotifsNotAdjacent(self):
        """
        The find method must find repeating identical motifs.
        """
        read = AARead('id', 'EKENLGYAAAAAAAAEGKENGY')
        finder = EukaryoticLinearMotif()
        result = list(finder.find(read))
        self.assertEqual([Landmark('EukaryoticLinearMotif', 'ELM', 0, 5,
                                   'DEG_APCC_KENBOX_2'),
                          Landmark('EukaryoticLinearMotif', 'ELM', 16, 5,
                                   'DEG_APCC_KENBOX_2')], result)

    def testFindTwoIdenticalMotifsAdjacent(self):
        """
        The find method must find two identical motifs next to each other.
        """
        read = AARead('id', 'EKENLGKENG')
        finder = EukaryoticLinearMotif()
        result = list(finder.find(read))
        self.assertEqual([Landmark('EukaryoticLinearMotif', 'ELM', 0, 5,
                                   'DEG_APCC_KENBOX_2'),
                          Landmark('EukaryoticLinearMotif', 'ELM', 5, 5,
                                   'DEG_APCC_KENBOX_2')], result)

    def testThreeDifferentMotifs(self):
        """
        The find method must find two different motifs.
        """
        read = AARead('id', 'EKENLGKENGWWWPWWP')
        finder = EukaryoticLinearMotif()
        result = list(finder.find(read))
        self.assertEqual([Landmark('EukaryoticLinearMotif', 'ELM', 0, 5,
                                   'DEG_APCC_KENBOX_2'),
                          Landmark('EukaryoticLinearMotif', 'ELM', 5, 5,
                                   'DEG_APCC_KENBOX_2'),
                          Landmark('EukaryoticLinearMotif', 'ELM', 10, 4,
                                   'LIG_CSL_BTD_1'),
                          Landmark('EukaryoticLinearMotif', 'ELM', 10, 7,
                                   'LIG_SH3_3')], result)

    def testMatchAlteration(self):
        """
        The find method must work on a regex that contains alternation (e.g.,
        [ABC]).
        """
        # Make sure we know what pattern is actually being matched.
        self.assertEqual('.P[TS]AP.', _DATABASE[22]['pattern'].pattern)
        self.assertEqual('LIG_PTAP_UEV_1', _DATABASE[22]['identifier'])

        finder = EukaryoticLinearMotif()

        result = list(finder.find(AARead('id', 'SPTAPS')))
        self.assertEqual([Landmark('EukaryoticLinearMotif', 'ELM', 0, 6,
                                   'LIG_PTAP_UEV_1')], result)

        result = list(finder.find(AARead('id', 'SPSAPS')))
        self.assertEqual([Landmark('EukaryoticLinearMotif', 'ELM', 0, 6,
                                   'LIG_PTAP_UEV_1')], result)

    def testLoadDatabase(self):
        """
        The EukaryoticLinearMotif database loading function must work and
        produce the expected result.
        """
        db = _loadDatabase()
        # The 49 is based on EukaryoticLinearMotif database elm-160713 of
        # 13-Jul-2016
        self.assertEqual(49, len(db))
        expectedKeys = {'pattern', 'identifier'}
        for motif in db:
            self.assertEqual(expectedKeys, set(motif.keys()))
