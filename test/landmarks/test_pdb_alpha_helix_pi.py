from unittest import TestCase

from dark.fasta_ss import SSAARead

from light.features import Landmark
from light.landmarks.pdb_alpha_helix_pi import PDB_AlphaHelix_pi


class TestPDB_AlphaHelix_pi(TestCase):
    """
    Tests for the PDB_AlphaHelix_pi.find function.
    """
    def testEmptyRead(self):
        """
        An empty read must not return any landmarks.
        """
        read = SSAARead('id', '', '')
        landmark = PDB_AlphaHelix_pi()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testNoLandmarks(self):
        """
        A read that has no landmarks must result in an empty list.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHHHHHHHH')
        landmark = PDB_AlphaHelix_pi()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testOneLandmarkAtStart(self):
        """
        If the structure sequence has the landmark at its beginning, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'IIIHHHHHH')
        landmark = PDB_AlphaHelix_pi()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix_pi', 'PDB-Api', 0, 3)],
                         list(result))

    def testOneLandmarkInMiddle(self):
        """
        If the structure sequence has the landmark in the middle, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHIIIHHHH')
        landmark = PDB_AlphaHelix_pi()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix_pi', 'PDB-Api', 2, 3)],
                         list(result))

    def testOneLandmarkAtEnd(self):
        """
        If the structure sequence has the landmark at its end, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHHHHHIII')
        landmark = PDB_AlphaHelix_pi()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix_pi', 'PDB-Api', 6, 3)],
                         list(result))

    def testTwoLandmarks(self):
        """
        The right landmarks must be returned when a structure sequence with
        two landmarks is given.
        """
        read = SSAARead(
            'id',
            'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            '-IIIIIIIIIIIIIIIIIIIIIIIIIIT-TT-TTSS---TTT-TTIIII-SS---')
        landmark = PDB_AlphaHelix_pi()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB AlphaHelix_pi', 'PDB-Api', 1, 26),
                          Landmark('PDB AlphaHelix_pi', 'PDB-Api', 45, 4)],
                         list(result))
