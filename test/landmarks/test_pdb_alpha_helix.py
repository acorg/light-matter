from unittest import TestCase

from dark.fasta_ss import SSAARead

from light.features import Landmark
from light.landmarks.pdb_alpha_helix import PDB_AlphaHelix


class TestPDB_AlphaHelix(TestCase):
    """
    Tests for the PDB_AlphaHelix.find function.
    """
    def testEmptyRead(self):
        """
        An empty read must not return any landmarks.
        """
        read = SSAARead('id', '', '')
        landmark = PDB_AlphaHelix()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testNoLandmarks(self):
        """
        A read that has no landmarks must result in an empty list.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'GGGGGGGGG')
        landmark = PDB_AlphaHelix()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testOneLandmarkAtStart(self):
        """
        If the structure sequence has the landmark at its beginning, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHHGGGGGG')
        landmark = PDB_AlphaHelix()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix', 'PDB-A', 0, 3)],
                         list(result))

    def testOneLandmarkInMiddle(self):
        """
        If the structure sequence has the landmark in the middle, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'GGHHHGGGG')
        landmark = PDB_AlphaHelix()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix', 'PDB-A', 2, 3)],
                         list(result))

    def testOneLandmarkAtEnd(self):
        """
        If the structure sequence has the landmark at its end, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'GGGGGGHHH')
        landmark = PDB_AlphaHelix()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix', 'PDB-A', 6, 3)],
                         list(result))

    def testTwoLandmarks(self):
        """
        The right landmarks must be returned when a structure sequence with
        two landmarks is given.
        """
        read = SSAARead(
            'id',
            'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            '-HHHHHHHHHHHHHHHHHHHHHHHHHHT-TT-TTSS---TTT-TTHHHH-SS---')
        landmark = PDB_AlphaHelix()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix', 'PDB-A', 1, 26),
                          Landmark('PDB AlphaHelix', 'PDB-A', 45, 4)],
                         list(result))
