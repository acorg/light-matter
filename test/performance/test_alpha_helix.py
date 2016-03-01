from unittest import TestCase

from dark.reads import AARead, Reads
from light.performance.alpha_helix import (
    analyzeAlphaHelix, analyzeAlphaHelices)
from light.performance.stats import Stats

# In the tests below, 'F' has "peak" hydrophobic, 'R' has "trough"
# hydropathy, and 'G' has an intermediate hydropathy value.


class TestAnalyzeAlphaHelix(TestCase):
    """
    Tests for the light.performance.alpha_helix.analyzeAlphaHelix function
    (not to be confused with analyzeAlphaHelices, below).
    """
    def testEmptyRead(self):
        """
        The analysis of an empty alpha helix must be as expected.
        """
        read = AARead('id', '')
        self.assertEqual(
            {
                'consecutivePeaks': [],
                'consecutiveTroughs': [],
                'extremaCount': 0,
                'final': None,
                'initial': None,
            },
            analyzeAlphaHelix(read))

    def testNoExtrema(self):
        """
        If an alpha helix has no extrema, the analysis must be as expected.
        """
        read = AARead('id', 'GGG')
        self.assertEqual(
            {
                'consecutivePeaks': [],
                'consecutiveTroughs': [],
                'extremaCount': 0,
                'final': None,
                'initial': None,
            },
            analyzeAlphaHelix(read))

    def testOnePeakContainingOneAA(self):
        """
        If an alpha helix has one peak made up of one amino acid, the analysis
        must be as expected.
        """
        read = AARead('id', 'GFG')
        self.assertEqual(
            {
                'consecutivePeaks': [1],
                'consecutiveTroughs': [],
                'extremaCount': 1,
                'final': 1,
                'initial': 1,
            },
            analyzeAlphaHelix(read))

    def testOnePeakContainingTwoAAs(self):
        """
        If an alpha helix has one peak made up of two amino acids, the analysis
        must be as expected.
        """
        read = AARead('id', 'GFFG')
        self.assertEqual(
            {
                'consecutivePeaks': [2],
                'consecutiveTroughs': [],
                'extremaCount': 2,
                'final': 1,
                'initial': 1,
            },
            analyzeAlphaHelix(read))

    def testOnePeakContainingTwoAAsWithIntermediate(self):
        """
        If an alpha helix has one peak made up of two amino acids separated by
        a non-peak AA, the analysis must be as expected.
        """
        read = AARead('id', 'GFGFGG')
        self.assertEqual(
            {
                'consecutivePeaks': [2],
                'consecutiveTroughs': [],
                'extremaCount': 2,
                'final': 2,
                'initial': 1,
            },
            analyzeAlphaHelix(read))

    def testTwoPeaksEachContainingOneAA(self):
        """
        If an alpha helix has two peaks, each made up of one amino acid, the
        analysis must be as expected.
        """
        read = AARead('id', 'GGFGRG')
        self.assertEqual(
            {
                'consecutivePeaks': [1],
                'consecutiveTroughs': [1],
                'extremaCount': 2,
                'final': 1,
                'initial': 2,
            },
            analyzeAlphaHelix(read))

    def testTwoPeaksEachContainingTwoAAs(self):
        """
        If an alpha helix has two peaks, each made up of two amino acids, the
        analysis must be as expected.
        """
        read = AARead('id', 'GGFFGRRG')
        self.assertEqual(
            {
                'consecutivePeaks': [2],
                'consecutiveTroughs': [2],
                'extremaCount': 4,
                'final': 1,
                'initial': 2,
            },
            analyzeAlphaHelix(read))

    def testTwoPeaksEachContainingTwoAAsWithIntermediate(self):
        """
        If an alpha helix has two peaks, each made up of two amino acids
        separated by a non-extrema AA, the analysis must be as expected.
        """
        read = AARead('id', 'GGFGFGRGRG')
        self.assertEqual(
            {
                'consecutivePeaks': [2],
                'consecutiveTroughs': [2],
                'extremaCount': 4,
                'final': 1,
                'initial': 2,
            },
            analyzeAlphaHelix(read))

    def testThreePeaksEachContainingOneAA(self):
        """
        If an alpha helix has three peaks, each made up of one amino acid, the
        analysis must be as expected.
        """
        read = AARead('id', 'GGGFGRFGGGG')
        self.assertEqual(
            {
                'consecutivePeaks': [1, 1],
                'consecutiveTroughs': [1],
                'extremaCount': 3,
                'final': 4,
                'initial': 3,
            },
            analyzeAlphaHelix(read))

    def testThreePeaksEachContainingMultipleAAsWithIntermediates(self):
        """
        If an alpha helix has three peaks, each made up of multiple amino
        acids with non-extrema intermediates, the analysis must be as expected.
        """
        read = AARead('id', 'GGGFGFFFGRRRGGRRFFGFGGGG')
        self.assertEqual(
            {
                'consecutivePeaks': [4, 3],
                'consecutiveTroughs': [5],
                'extremaCount': 12,
                'final': 4,
                'initial': 3,
            },
            analyzeAlphaHelix(read))

    def testZeroInitialAndFinal(self):
        """
        If an alpha helix has no leading non-extrema AAs, the 'initial'
        and 'final' values in the analysis must be zero.
        """
        read = AARead('id', 'F')
        self.assertEqual(
            {
                'consecutivePeaks': [1],
                'consecutiveTroughs': [],
                'extremaCount': 1,
                'final': 0,
                'initial': 0,
            },
            analyzeAlphaHelix(read))


class TestAnalyzeAlphaHelices(TestCase):
    """
    Tests for the light.performance.alpha_helix.analyzeAlphaHelices function
    (not to be confused with analyzeAlphaHelix, above).
    """
    def testResultDictValuesAreStatsInstances(self):
        """
        In the values analysis C{dict} returned by analyzeAlphaHelices, all
        values must be C{Stats} instances.
        """
        reads = Reads()
        for value in analyzeAlphaHelices(reads).values():
            self.assertIsInstance(value, Stats)

    def testNoHelices(self):
        """
        The analysis of an empty set of alpha helices must be as expected.
        """
        reads = Reads()
        self.assertEqual(
            {
                'length': [],
                'initial': [],
                'final': [],
                'extremaCount': [],
                'consecutivePeaks': [],
                'consecutiveTroughs': [],
                'noExtremaLength': [],
            },
            analyzeAlphaHelices(reads))

    def testOneHelix(self):
        """
        The analysis of a single alpha helix must be as expected.
        """
        reads = Reads()
        reads.add(AARead('id', 'F'))
        self.assertEqual(
            {
                'length': [1],
                'initial': [0],
                'final': [0],
                'extremaCount': [1],
                'consecutivePeaks': [1],
                'consecutiveTroughs': [],
                'noExtremaLength': [],
            },
            analyzeAlphaHelices(reads))

    def testTwoIdenticalHelices(self):
        """
        The analysis of two identical alpha helices must be as expected.
        """
        reads = Reads()
        reads.add(AARead('id', 'F'))
        reads.add(AARead('id', 'F'))
        self.assertEqual(
            {
                'length': [1, 1],
                'initial': [0, 0],
                'final': [0, 0],
                'extremaCount': [1, 1],
                'consecutivePeaks': [1, 1],
                'consecutiveTroughs': [],
                'noExtremaLength': [],
            },
            analyzeAlphaHelices(reads))

    def testTwoDifferingHelices(self):
        """
        The analysis of two differing alpha helices must be as expected.
        """
        reads = Reads()
        reads.add(AARead('id', 'FG'))
        reads.add(AARead('id', 'GRR'))
        self.assertEqual(
            {
                'length': [2, 3],
                'initial': [0, 1],
                'final': [1, 0],
                'extremaCount': [1, 2],
                'consecutivePeaks': [1],
                'consecutiveTroughs': [2],
                'noExtremaLength': [],
            },
            analyzeAlphaHelices(reads))

    def testNoExtrema(self):
        """
        The analysis of a helix with no extrema must be correct.
        """
        reads = Reads()
        reads.add(AARead('id', 'GG'))
        self.assertEqual(
            {
                'length': [2],
                'initial': [],
                'final': [],
                'extremaCount': [0],
                'consecutivePeaks': [],
                'consecutiveTroughs': [],
                'noExtremaLength': [2],
            },
            analyzeAlphaHelices(reads))

    def testFourDifferingHelices(self):
        """
        The analysis of four differing alpha helices must be as expected.
        """
        reads = Reads()
        reads.add(AARead('id1', 'GGFFFRRFFFG'))
        reads.add(AARead('id2', 'GRRFFFFGGGRRRGRGRGGFFGGG'))
        reads.add(AARead('id3', 'RRGRRRFFGFFGGGRRR'))
        reads.add(AARead('id4', 'GGG'))
        self.assertEqual(
            {
                'length': [11, 24, 17, 3],
                'initial': [2, 1, 0],
                'final': [1, 3, 0],
                'extremaCount': [8, 13, 12, 0],
                'consecutivePeaks': [3, 3,     # id1
                                     4, 2,     # id2
                                     4],       # id3
                'consecutiveTroughs': [2,      # id1
                                       2, 5,   # id2
                                       5, 3],  # id3
                'noExtremaLength': [3],        # id4
            },
            analyzeAlphaHelices(reads))
