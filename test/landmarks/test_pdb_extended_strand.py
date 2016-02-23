from unittest import TestCase

from dark.fasta_ss import SSAARead

from light.features import Landmark
from light.landmarks.pdb_extended_strand import PDB_ExtendedStrand


class TestPDB_ExtendedStrand(TestCase):
    """
    Tests for the PDB_ExtendedStrand.find function.
    """
    def testEmptyRead(self):
        """
        An empty read must not return any landmarks.
        """
        read = SSAARead('id', '', '')
        landmark = PDB_ExtendedStrand()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testNoLandmarks(self):
        """
        A read that has no landmarks must result in an empty list.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHHHHHHHH')
        landmark = PDB_ExtendedStrand()
        result = landmark.find(read)
        self.assertEqual([], list(result))

    def testOneLandmarkAtStart(self):
        """
        If the structure sequence has the landmark at its beginning, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'EEEHHHHHH')
        landmark = PDB_ExtendedStrand()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB ExtendedStrand', 'PDB-ES', 0, 3)],
                         list(result))

    def testOneLandmarkInMiddle(self):
        """
        If the structure sequence has the landmark in the middle, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHEEEHHHH')
        landmark = PDB_ExtendedStrand()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB ExtendedStrand', 'PDB-ES', 2, 3)],
                         list(result))

    def testOneLandmarkAtEnd(self):
        """
        If the structure sequence has the landmark at its end, the right
        landmark must be returned.
        """
        read = SSAARead('id', 'RFRFRFRFR', 'HHHHHHEEE')
        landmark = PDB_ExtendedStrand()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB ExtendedStrand', 'PDB-ES', 6, 3)],
                         list(result))

    def testTwoLandmarks(self):
        """
        The right landmarks must be returned when a structure sequence with
        two landmarks is given.
        """
        read = SSAARead(
            'id',
            'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            '-EEEEEEEEEEEEEEEEEEEEEEEEEET-TT-TTSS---TTT-TTEEEE-SS---')
        landmark = PDB_ExtendedStrand()
        result = landmark.find(read)
        self.assertEqual([Landmark('PDB ExtendedStrand', 'PDB-ES', 1, 26),
                          Landmark('PDB ExtendedStrand', 'PDB-ES', 45, 4)],
                         list(result))
