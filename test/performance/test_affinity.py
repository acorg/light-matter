from unittest import TestCase
import numpy as np

from dark.reads import Reads, AARead

from light.parameters import FindParameters
from light.performance.affinity import affinityMatrix
from light.landmarks import ALL_LANDMARK_CLASSES
from light.trig import ALL_TRIG_CLASSES


class TestAffinityMatrix(TestCase):
    """
    Tests for the light.performance.affinity.affinityMatrix function.
    """

    def testNoReads(self):
        """
        If affinityMatrix is called with no reads and no subjects, an empty
        matrix must be returned.
        """
        reads = Reads()
        matrix = affinityMatrix(reads, landmarks=['AlphaHelix'])
        self.assertTrue(np.array_equal(np.zeros((0, 0)), matrix))

    def testSequenceWithNoFeaturesAgainstItself(self):
        """
        If affinityMatrix is called with a read that is also the only subject
        in the database, and the read has no features, a matrix with just a
        single 0.0 value must be returned.
        """
        reads = Reads()
        read = AARead('id1', 'AAA')
        reads.add(read)
        matrix = affinityMatrix(reads, landmarks=['AlphaHelix'],
                                computeDiagonal=True)
        self.assertEqual([[0.0]], matrix)

    def testOneSequenceSpecificDiagonalValue(self):
        """
        If affinityMatrix is called with a single read and a specific
        diagonal value, that diagonal value must be in the result.
        """
        reads = Reads()
        read = AARead('id1', 'AAA')
        reads.add(read)
        matrix = affinityMatrix(reads, landmarks=['AlphaHelix'],
                                diagonalValue=2.0)
        self.assertEqual([[2.0]], matrix)

    def testSequenceWithFeaturesAgainstItself(self):
        """
        If affinityMatrix is called with a read that is also the only subject
        in the database, and the read has features, a matrix with just a
        single 1.0 value must be returned.
        """
        reads = Reads()
        read = AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF')
        reads.add(read)
        matrix = affinityMatrix(reads, landmarks=['AlphaHelix'],
                                computeDiagonal=True)
        self.assertEqual([[1.0]], matrix)

    def testZeroByThree(self):
        """
        If affinityMatrix is called with no reads and three subjects, the
        resulting matrix must be 0x3.
        """
        reads = Reads()
        subjects = Reads()
        subjects.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id3', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id4', 'FRRRFRRRFAAAFRRRFRRRF'))
        matrix = affinityMatrix(reads, landmarks=['AlphaHelix'],
                                subjects=subjects, computeDiagonal=True)
        self.assertTrue(np.array_equal(np.zeros((0, 3)), matrix))

    def testOneByZero(self):
        """
        If affinityMatrix is called with no reads and three subjects, the
        resulting matrix must be 1x0.
        """
        reads = Reads()
        read = AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF')
        reads.add(read)
        subjects = Reads()
        matrix = affinityMatrix(reads, landmarks=['AlphaHelix'],
                                subjects=subjects, computeDiagonal=True)
        self.assertTrue(np.array_equal(np.zeros((1, 0)), matrix))

    def testOneByThree(self):
        """
        If affinityMatrix is called with a read and the database has three
        subjects, the resulting matrix must be 1x3.
        """
        reads = Reads()
        read = AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF')
        reads.add(read)
        subjects = Reads()
        subjects.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id3', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id4', 'FRRRFRRRFAAAFRRRFRRRF'))
        matrix = affinityMatrix(reads, landmarks=['AlphaHelix'],
                                subjects=subjects, computeDiagonal=True)
        self.assertTrue(np.array_equal([[1.0, 1.0, 1.0]], matrix))

    def testTwoByTwoWithProgressFunction(self):
        """
        If affinityMatrix is called with two reads and the database has two
        subjects, and a progress function is passed, the progress function
        must be called as expected.
        """
        reads = Reads()
        reads.add(AARead('id1', 'FRRRFRRRFAAAFRRRFRRRF'))
        reads.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects = Reads()
        subjects.add(AARead('id3', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id4', 'FRRRFRRRFAAAFRRRFRRRF'))

        output = []

        def progress(i, query):
            output.append((i, query.id))

        affinityMatrix(reads, landmarks=['AlphaHelix'],
                       subjects=subjects, computeDiagonal=True,
                       progressFunc=progress)

        self.assertEqual([(0, 'id1'), (1, 'id2')], output)

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
        subjects.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id3', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id4', 'FRRRFRRRFAAAFRRRFRRRF'))
        matrix = affinityMatrix(reads, landmarks=['AlphaHelix'],
                                subjects=subjects, computeDiagonal=True)
        self.assertTrue(np.array_equal(
            [
                [1.0, 1.0, 1.0],
                [1.0, 1.0, 1.0]
            ],
            matrix)
        )

    def _checkSymmetry(self, sequences, findParams, symmetric=False, **kwargs):
        """
        Create an affinity matrix for a set of sequences and check its
        symmetry.

        @param sequences: A C{list} of C{AARead} instances.
        @param findParams: A {light.parameters.FindParameters} instance.
        @param symmetric: If C{True}, pass symmetric=True to the affinityMatrix
            function, allowing it to speed up the calculation by assuming
            scores are symmetric. We still check that the result actually is
            symmetric.
        @param kwargs: See
            C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
            additional keywords, all of which are optional.
        """
        matrix = affinityMatrix(sequences, findParams, symmetric=symmetric,
                                **kwargs)

        for i in range(len(sequences)):

            # Test the diagonal score of each sequence against itself is 1.0.
            self.assertEqual(
                1.0, matrix[i][i],
                'Diagonal entry (%d, %d) for %s against itself has non-1.0 '
                'score of %f.' % (i, i, sequences[i].id, matrix[i][i]))

            # Test that off-diagonal score pairs are identical.
            for j in range(i + 1, len(sequences)):
                self.assertEqual(
                    matrix[i][j], matrix[j][i],
                    'Off-diagonal entries (%d, %d) and (%d, %d) for %s '
                    'against %s have unequal scores %f and %f.' %
                    (i, j, j, i, sequences[i].id, sequences[j].id,
                     matrix[i][j], matrix[j][i]))

    def testSandraSymmetry_235(self):
        """
        Make sure we get a symmetric affinity matrix on a few of the sequences
        received from Sandra Junglen on March 13, 2015.

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

        findParams = FindParameters(significanceFraction=0.05)
        self._checkSymmetry(
            sequences, findParams, distanceBase=1.0,
            landmarks=ALL_LANDMARK_CLASSES,
            trigPoints=ALL_TRIG_CLASSES,
            limitPerLandmark=50, minDistance=1, maxDistance=100,
            symmetric=False)

    def testSandraSymmetryWithSymmetricSpeedup_235(self):
        """
        Make sure we get a symmetric affinity matrix on a few of the sequences
        received from Sandra Junglen on March 13, 2015 if we pass no value for
        'symmetric' (which defaults to True) to affinityMatrix.

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

        findParams = FindParameters(significanceFraction=0.05)
        self._checkSymmetry(
            sequences, findParams, distanceBase=1.0,
            landmarks=ALL_LANDMARK_CLASSES,
            trigPoints=ALL_TRIG_CLASSES,
            limitPerLandmark=50, minDistance=1, maxDistance=100)

    def testSandraSymmetry_259(self):
        """
        Make sure we get a symmetric affinity matrix on two of the sequences
        received from Sandra Junglen on March 13, 2015.

        The sequences below are the ones that caused the non-symmetric
        scores issue in https://github.com/acorg/light-matter/issues/259
        """
        sequences = [
            # Read index 8.
            AARead('RGSV',
                   ('HVCIFRKNQHGGLREIYVLNIYERIVQKCVEDLARAILSVVPSETMTHPKNKF'
                    'QIPNKHNIAARKEFGDSYFTVCTSDDASKWNQGHHVSKFITILVRILPKFWHG'
                    'FIVRALQLWFHKRLFLGDDLLRLFCANDVLNTTDEKVKKVHEVFKGREVAPWM'
                    'TRGMTYIETESGFMQGILHYISSLFHAIFLEDLAERQKKQLPQMARIIQPDNE'
                    'SNVIIDCMESSDDSSMMISFSTKSMNDRQTFAMLLLVDRAFSLKEYYGDMLGI'
                    'YKSIKSTTGTIFMMEFNIEFFFAGDTHRPTIRWVNAALN')),

            # Read index 47.
            AARead('InfluenzaC',
                   ('TINTMAKDGERGKLQRRAIATPGMIVRPFSKIVETVAQKICEKLKESGLPVGG'
                    'NEKKAKLKTTVTSLNARMNSDQFAVNITGDNSKWNECQQPEAYLALLAYITKD'
                    'SSDLMKDLCSVAPVLFCNKFVKLGQGIRLSNKRKTKEVIIKAEKMGKYKNLMR'
                    'EEYKNLFEPLEKYIQKDVCFLPGGMLMGMFNMLSTVLGVSTLCYMDEELKAKG'
                    'CFWTGLQSSDDFVLFAVASNWSNIHWTIRRFNAVCKLIGINMSLEKSYGSLPE'
                    'LFEFTSMFFDGEFVSNLAMELPAFT')),
        ]

        findParams = FindParameters(significanceFraction=0.01)
        self._checkSymmetry(
            sequences, findParams, distanceBase=1.025,
            landmarks=['GOR4AlphaHelix', 'GOR4Coil'],
            trigPoints=['Peaks', 'Troughs'],
            limitPerLandmark=50, minDistance=1, maxDistance=100,
            symmetric=False)
