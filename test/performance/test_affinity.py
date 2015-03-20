from unittest import TestCase

from dark.reads import Reads, AARead

from light.performance.affinity import affinityMatrix
from light.landmarks import ALL_LANDMARK_CLASSES
from light.trig import ALL_TRIG_CLASSES


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

    def testSandraSymmetry(self):
        """
        Make sure we get a symmetric affinity matrix on a few of the sequences
        received from Sandra Junglen on March 13, 2015. We need 1.0 scores on
        the diagonal and pairs of corresponding off-diagonal scores must be
        equal.

        The sequences below are the ones that caused the non-symmetric
        scores issue in https://github.com/acorg/light-matter/issues/235
        """
        sequences = [
            # Read index 3.
            AARead('BUNV', ('SFTFFNKGQKTAKDREIFVGEFEAKMCMYVVERISKERCKLNTDE'
                            'MISEPGDSKLKILEKKAEEEIRYIVERTKDSIIKGDPSKALKLEI'
                            'NADMSKWSAQDVFYKYFWLIAMDPILYPAEKTRILYFMCNYMQKL'
                            'LILPDDLIANILDQKRPYNDDLILEMTNGLNYNYVQIKRNWLQGN'
                            'FNYISSYVHSCAMLVYKDILKECMKLLDGDCLINSMVHSDDNQTS'
                            'LAIIQNKVSDQIVIQYAANTFESVCLTFGCQANMKKTYITHTCKE'
                            'FVSLFNLHGEPLSVFGRFLLPSVG')),

            # Read index 24.
            AARead('LACV', ('YFTFFNKGQKTSKDREIFVGEYEAKMCMYAVERIAKERCKLNPDE'
                            'MISEPGDGKLKVLEQKSEQEIRFLVETTRQKNREIDEAIEALAAE'
                            'GYESNLEKIEKLSLGKAKGLKMEINADMSKWSAQDVFYKYFWLIA'
                            'LDPILYPQEKERILYFMCNYMDKELILPDELLFNLLDQKVAYQND'
                            'IIATMTNQLNSNTVLIKRNWLQGNFNYTSSYVHSCAMSVYKEILK'
                            'EAITLLDGSILVNSLVHSDDNQTSITIVQDKMENDKIIDFAMKEF'
                            'ERACLTFGCQANMKKTYVTNCIKEFVSLFNLYGEPFSIYGRFLLT'
                            'SVG')),

            # Read index 48.
            AARead('WYOV', ('TFTFFNKGQKTAKDREIFVGEFEAKMCMYVVERIAKERCKLNSDE'
                            'MISEPGDAKLKILEQKAEQELRFIVERTKDKFLKGDPCKALKMEI'
                            'NADMSKWSAQDVFYKYFWLIAMDPILYPKEKYRILFFMCNYLQKV'
                            'LVLPDELIGNILDQKKTYNNDIILEGTDFLHQNYVNIRRNWLQGN'
                            'FNYLSSYIHTCAMSVFKDILKEVSYLLDGDVLVNSMVHSDDNQTS'
                            'ITYVQNKIEESVLINHGLKTFETVCLTFGCQANMKKTYLTHNIKE'
                            'FVSLFNIHGEPMSVYGRFLLPSVG')),
        ]

        matrix = affinityMatrix(
            sequences, significanceFraction=0.05, distanceBase=1.0,
            landmarkNames=[cls.NAME for cls in ALL_LANDMARK_CLASSES],
            trigPointNames=[cls.NAME for cls in ALL_TRIG_CLASSES],
            limitPerLandmark=50, minDistance=1, maxDistance=100,
            subjects=sequences)

        for i in xrange(len(sequences)):

            # Test the diagonal score of a sequence against itself is one.
            self.assertEqual(
                1.0, matrix[i][i],
                'Diagonal entry (%d, %d) for %s against itself has non-1.0 '
                'score of %f.' % (i, i, sequences[i].id, matrix[i][i]))

            # Test that off-diagonal score pairs are identical.
            for j in xrange(i + 1, len(sequences)):
                self.assertEqual(
                    matrix[i][j], matrix[j][i],
                    'Off-diagonal entries (%d, %d) and (%d, %d) for %s '
                    'against %s have unequal scores %f and %f.' %
                    (i, j, j, i, sequences[i].id, sequences[j].id,
                     matrix[i][j], matrix[j][i]))
