from unittest import TestCase

from dark.reads import Reads, AARead

from light.performance.affinity import affinityMatrix


class TestAffinityMatrix(TestCase):
    """
    Tests for the light.performance.affinity.affinityMatrix function.
    """

    def testNoReads(self):
        """
        If affinityMatrix is called with no reads, an empty matrix must be
        returned.
        """
        reads = Reads()
        matrix = affinityMatrix(reads, landmarkNames=['AlphaHelix'])
        self.assertEqual([], matrix)

    def testEmptyDatabase(self):
        """
        If affinityMatrix is called with reads but the database is empty, an
        empty matrix must be returned.
        """
        reads = Reads()
        matrix = affinityMatrix(reads, landmarkNames=['AlphaHelix'],
                                subjects=reads)
        self.assertEqual([], matrix)

    def testSequenceWithNoFeaturesAgainstItself(self):
        """
        If affinityMatrix is called with a read that is also the only subject
        in the database, and the read has no features, a matrix with just a
        single 0.0 value must be returned.
        """
        reads = Reads()
        read = AARead('id1', 'AAA')
        reads.add(read)
        matrix = affinityMatrix(reads, landmarkNames=['AlphaHelix'],
                                subjects=reads)
        self.assertEqual([[0.0]], matrix)

    def testSequenceWithFeaturesAgainstItself(self):
        """
        If affinityMatrix is called with a read that is also the only subject
        in the database, and the read has features, a matrix with just a
        single 1.0 value must be returned.
        """
        reads = Reads()
        read = AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF')
        reads.add(read)
        matrix = affinityMatrix(reads, landmarkNames=['AlphaHelix'],
                                subjects=reads)
        self.assertEqual([[1.0]], matrix)

    def testOneByThree(self):
        """
        If affinityMatrix is called with a read and the database has three
        subjects, the resulting matrix must be 1x3.
        """
        reads = Reads()
        read = AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF')
        reads.add(read)
        subjects = Reads()
        subjects.add(read)
        subjects.add(read)
        subjects.add(read)
        matrix = affinityMatrix(reads, landmarkNames=['AlphaHelix'],
                                subjects=subjects)
        self.assertEqual([[1.0, 1.0, 1.0]], matrix)

    def testTwoByThree(self):
        """
        If affinityMatrix is called with two reads and the database has three
        subjects, the resulting matrix must be 2x3.
        """
        reads = Reads()
        read = AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF')
        reads.add(read)
        reads.add(read)
        subjects = Reads()
        subjects.add(read)
        subjects.add(read)
        subjects.add(read)
        matrix = affinityMatrix(reads, landmarkNames=['AlphaHelix'],
                                subjects=subjects)
        self.assertEqual([
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0]
        ], matrix)
