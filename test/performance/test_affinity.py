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
        subjects.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id3', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id4', 'FRRRFRRRFAAAFRRRFRRRF'))
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
        subjects.add(AARead('id2', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id3', 'FRRRFRRRFAAAFRRRFRRRF'))
        subjects.add(AARead('id4', 'FRRRFRRRFAAAFRRRFRRRF'))
        matrix = affinityMatrix(reads, landmarkNames=['AlphaHelix'],
                                subjects=subjects)
        self.assertEqual([
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0]
        ], matrix)

    def _checkSymmetry(self, sequences, **kwargs):
        """
        Create an affinity matrix for a set of sequences and check its
        symmetry.

        @param sequences: A C{list} of C{AARead} instances.
        @param kwargs: See
            C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
            additional keywords, all of which are optional.
        """
        matrix = affinityMatrix(sequences, subjects=sequences, **kwargs)

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

        self._checkSymmetry(
            sequences, significanceFraction=0.05, distanceBase=1.0,
            landmarkNames=[cls.NAME for cls in ALL_LANDMARK_CLASSES],
            trigPointNames=[cls.NAME for cls in ALL_TRIG_CLASSES],
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

        self._checkSymmetry(
            sequences, significanceFraction=0.01, distanceBase=1.025,
            landmarkNames=['GOR4AlphaHelix', 'GOR4Coil'],
            trigPointNames=['Peaks', 'Troughs'],
            limitPerLandmark=50, minDistance=1, maxDistance=100)
