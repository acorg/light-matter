from unittest import TestCase

from dark.fasta_ss import SSAARead

from light.features import Landmark
from light.landmarks.pdb_alpha_helix_combined import (
    PDB_AlphaHelix_combined, combineHelices)


class TestPDB_AlphaHelix_combined(TestCase):
    """
    Tests for the PDB_AlphaHelix_combined.find function.
    """
    def testEmptyRead(self):
        """
        An empty read must not return any landmarks.
        """
        read = SSAARead('id', '', '')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testNoLandmarks(self):
        """
        A read that has no landmarks must result in an empty list.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'SSSSSSSSS')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testAlphaHelixMustBeFound(self):
        """
        If the stucture sequence contains alpha helix ('H'), the right landmark
        must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHHSSSSSS')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual([
                         Landmark('PDB AlphaHelix_combined', 'PDB-AC', 0, 3)
                         ],
                         list(result))

    def testAlphaHelix310MustBeFound(self):
        """
        If the stucture sequence contains alpha helix 3 10 ('G'), the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'GGGSSSSSS')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual([
                         Landmark('PDB AlphaHelix_combined', 'PDB-AC', 0, 3)
                         ],
                         list(result))

    def testAlphaHelixPiMustBeFound(self):
        """
        If the stucture sequence contains alpha helix pi ('I'), the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'IIISSSSSS')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual([
                         Landmark('PDB AlphaHelix_combined', 'PDB-AC', 0, 3)
                         ],
                         list(result))

    def testAlphaHelixMixtureMustBeFound(self):
        """
        If the stucture sequence contains all alpha helix types ('H', 'G',
        'I'), the right landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HGISSSSSS')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual(
            [Landmark('PDB AlphaHelix_combined', 'PDB-AC', 0, 3)],
            list(result))

    def testOneLandmarkAtStart(self):
        """
        If the structure sequence has the landmark at its beginning, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHHSSSSSS')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual(
            [Landmark('PDB AlphaHelix_combined', 'PDB-AC', 0, 3)],
            list(result))

    def testOneLandmarkInMiddle(self):
        """
        If the structure sequence has the landmark in the middle, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'SSHHHSSSS')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual(
            [Landmark('PDB AlphaHelix_combined', 'PDB-AC', 2, 3)],
            list(result))

    def testOneLandmarkAtEnd(self):
        """
        If the structure sequence has the landmark at its end, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'SSSSSSHHH')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual(
            [Landmark('PDB AlphaHelix_combined', 'PDB-AC', 6, 3)],
            list(result))

    def testTwoLandmarks(self):
        """
        The right landmarks must be returned when a structure sequence with
        two landmarks is given.
        """
        read = SSAARead(
            'id',
            'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            '-HHHHHHHHHHHHHHHGGHHHHHIIHHT-TT-TTSS---TTT-TTHHHH-SS---')
        landmark = PDB_AlphaHelix_combined()
        result = landmark.find(read)
        self.assertEqual(
            [
                Landmark('PDB AlphaHelix_combined', 'PDB-AC', 1, 26),
                Landmark('PDB AlphaHelix_combined', 'PDB-AC', 45, 4),
            ], list(result))


class TestCombineHelices(TestCase):
    """
    Tests for the combineHelices function.
    """
    def testEmptyStructure(self):
        """
        If an empty structure sequence is passed in, the correct result must be
        returned.
        """
        result = combineHelices('')
        self.assertEqual('', result)

    def testBeginning(self):
        """
        If a helix is at the beginning of a sequence, the correct result must
        be returned.
        """
        result = combineHelices('HHHSSSSS')
        self.assertEqual('KKKCCCCC', result)

    def testMiddle(self):
        """
        If a helix is in the middle of a sequence, the correct result must be
        returned.
        """
        result = combineHelices('SSHHHSS')
        self.assertEqual('CCKKKCC', result)

    def testEnd(self):
        """
        If a helix is at the end of a sequence, the correct result must be
        returned.
        """
        result = combineHelices('SSSSHHH')
        self.assertEqual('CCCCKKK', result)

    def testAlphaHelix(self):
        """
        If a string of AlphaHelix ('H') is given, the correct string must be
        returned.
        """
        result = combineHelices('HHHHH')
        self.assertEqual('KKKKK', result)

    def testAlphaHelix310(self):
        """
        If a string of AlphaHelix_3_10 ('G') is given, the correct string must
        be returned.
        """
        result = combineHelices('GGGGG')
        self.assertEqual('KKKKK', result)

    def testAlphaHelixPi(self):
        """
        If a string of AlphaHelix_pi ('I') is given, the correct string must be
        returned.
        """
        result = combineHelices('IIIII')
        self.assertEqual('KKKKK', result)

    def testAllTypes(self):
        """
        If a string of a mixture of all alpha helix structure types ('H', 'G',
        'I') is given, the correct string must be returned.
        """
        result = combineHelices('HGIGH')
        self.assertEqual('KKKKK', result)
