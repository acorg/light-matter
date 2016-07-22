from unittest import TestCase

from dark.reads import SSAARead

from light.performance.overlap import CalculateOverlap


class TestCalculateOverlap(TestCase):
    """
    Tests for the CalculateOverlap class.
    """
    def testGetFeatures(self):
        """
        Calling getFeatures on an SSAARead must give the correct result.
        """
        sequence = 'SMEQVAMELRLTELTRLLRSVLDQLQDKDPARIFAQPVSLKEVPDYLDHIKHPMD'
        structure = '-HHHHHHHHHHHHHHHHHHHHHHHHHHT-TT-TTSS---TTT-TTHHHH-SS---'
        ssAARead = SSAARead('5AMF', sequence, structure)
        co = CalculateOverlap()
        features, intersection, union = co.getFeatures(ssAARead)
        expected = {
            'AC AlphaHelix_pi': set(),
            'AlphaHelix': set(),
            'BetaTurn': set(),
            'AC AlphaHelix_3_10': set(),
            'AminoAcidsLm': set(),
            'GOR4Coil': {
                0, 1, 2, 3, 4, 5, 6, 25, 26, 27, 28, 29, 34, 35, 36, 37, 38,
                39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
            },
            'IndividualPeaks': {
                9, 31,
            },
            'PDB AlphaHelix_combined': {
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                19, 20, 21, 22, 23, 24, 25, 26, 45, 46, 47, 48
            },
            'IndividualTroughs': {
                49, 10,
            },
            'ClusterAlphaHelix': {
                33, 34, 35, 36, 10, 11, 12, 13, 14, 15, 16, 17,
            },
            'Troughs': {
                1, 33, 4, 37, 6, 39, 8, 10, 42, 13, 46, 49, 21, 53, 24, 30,
            },
            'Peaks': {
                3, 35, 5, 38, 7, 40, 9, 12, 44, 15, 47, 18, 50, 22, 27, 31,
            },
            'AminoAcids': set(),
            'AC AlphaHelix': {
                32, 33, 34, 35, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 31,
            },
            'AC AlphaHelix_combined': {
                32, 33, 34, 35, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 31,
            },
            'GOR4BetaStrand': set(),
            'EukaryoticLinearMotif': {
                9, 10, 11, 12, 13, 14, 15, 16, 17, 35, 36, 37, 38, 39, 40, 41,
                45, 46, 47, 48,
            },
            'AC ExtendedStrand': {
                18, 19, 20, 21, 22,
            },
            'PDB AlphaHelix_3_10': set(),
            'BetaStrand': set(),
            'THAlphaHelix': {
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
                35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                51, 52, 53, 54,
            },
            'GOR4AlphaHelix': {
                7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                23, 24, 30, 31, 32, 33,
            },
            'AlphaHelix_pi': set(),
            'PDB AlphaHelix': {
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                19, 20, 21, 22, 23, 24, 25, 26, 45, 46, 47, 48,
            },
            'Prosite': {
                32, 38, 39, 40, 41, 14, 15, 16, 19, 20, 21, 22, 30, 31,
            },
            'AlphaHelix_3_10': {
                21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
            },
            'PDB AlphaHelix_pi': set(),
            'PDB ExtendedStrand': set(),
        }

        self.assertEqual(expected, features)

        self.assertEqual(
            {
                49, 10,
            },
            intersection[frozenset(('IndividualTroughs', 'Troughs'))])

        self.assertEqual(
            {
                1, 33, 4, 37, 6, 39, 8, 10, 42, 13, 46, 49, 21, 53, 24, 30,
            },
            union[frozenset(('IndividualTroughs', 'Troughs'))])

        # Note that the following don't test much. There are 28 features
        # examined by default by CalculateOverlap. So there are 28 * 27 / 2
        # = 378 pairs of features. So these two tests are just testing that
        # all pairs of features are present in the returned dicts.
        self.assertEqual(378, len(intersection))
        self.assertEqual(378, len(union))
