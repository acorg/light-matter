from unittest import TestCase

from dark.fasta_ss import SSAARead

from light.features import Landmark
from light.landmarks.pdb_alpha_helix_3_10 import PDB_AlphaHelix_3_10


class TestPDB_AlphaHelix_3_10(TestCase):
    """
    Tests for the PDB_AlphaHelix_3_10.find function.
    """
    def testEmptyRead(self):
        """
        An empty read must not return any landmarks.
        """
        read = SSAARead('id', '', '')
        landmark = PDB_AlphaHelix_3_10()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testNoLandmarks(self):
        """
        A read that has no landmarks must result in an empty list.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHHHHHHHH')
        landmark = PDB_AlphaHelix_3_10()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testOneLandmarkAtStart(self):
        """
        If the structure sequence has the landmark at its beginning, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'GGGHHHHHH')
        landmark = PDB_AlphaHelix_3_10()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix_3_10', 'PDB-A310', 0, 3)],
                         list(result))

    def testOneLandmarkInMiddle(self):
        """
        If the structure sequence has the landmark in the middle, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHGGGHHHH')
        landmark = PDB_AlphaHelix_3_10()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix_3_10', 'PDB-A310', 2, 3)],
                         list(result))

    def testOneLandmarkAtEnd(self):
        """
        If the structure sequence has the landmark at its end, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHHHHHGGG')
        landmark = PDB_AlphaHelix_3_10()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix_3_10', 'PDB-A310', 6, 3)],
                         list(result))

    def testTwoLandmarks(self):
        """
        The right landmarks must be returned when a structure sequence with
        two landmarks is given.
        """
        read = SSAARead(
            'id',
            'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            '-GGGGGGGGGGGGGGGGGGGGGGGGGGT-TT-TTSS---TTT-TTGGGG-SS---')
        landmark = PDB_AlphaHelix_3_10()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix_3_10', 'PDB-A310', 1, 26),
                          Landmark('PDB AlphaHelix_3_10', 'PDB-A310', 45, 4)],
                         list(result))
