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

        self.assertEqual(
            {
                'AlphaHelix': set(),
                'AlphaHelix_3_10': {
                    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                },
                'AlphaHelix_pi': set(),
                'AminoAcidsLm': set(),
                'BetaStrand': set(),
                'BetaTurn': set(),
                'GOR4AlphaHelix': {
                    5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                    21, 22, 23, 24, 29, 30, 31, 32, 33,
                },
                'GOR4BetaStrand': {
                    53, 54,
                },
                'GOR4Coil': {
                    0, 1, 2, 3, 4, 25, 26, 27, 28, 34, 35, 36, 37, 38, 39, 40,
                    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52,
                },
                'Prosite': {
                    32, 38, 39, 40, 41, 14, 15, 16, 19, 20, 21, 22, 30, 31,
                },
                'AminoAcids': set(),
                'IndividualPeaks': {
                    9, 31,
                },
                'IndividualTroughs': {
                    49, 10,
                },
                'Peaks': {
                    3, 35, 5, 38, 7, 40, 9, 12, 44, 15, 47, 18, 50, 22, 27, 31,
                },
                'Troughs': {
                    1, 33, 4, 37, 6, 39, 8, 10, 42, 13, 46, 49, 21, 53, 24, 30,
                },
                '-': {
                    0, 36, 37, 38, 42, 49, 52, 53, 54, 28, 31,
                },
                'E': set(),
                'G': set(),
                'H': {
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                    18, 19, 20, 21, 22, 23, 24, 25, 26, 45, 46, 47, 48,
                },
                'I': set(),
                'S': {
                    51, 34, 35, 50,
                },
                'T': {
                    32, 33, 39, 40, 41, 43, 44, 27, 29, 30,
                }
            },
            features)

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

        # Note that the following don't test much. There are 19 features
        # examined by default by CalculateOverlap. So there are 19 * 18 / 2
        # = 171 pairs of features. So these two tests are just testing that
        # all pairs of features are present in the returned dicts.
        self.assertEqual(171, len(intersection))
        self.assertEqual(171, len(union))
