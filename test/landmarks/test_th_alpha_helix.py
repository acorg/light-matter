from unittest import TestCase

from dark.aa import PROPERTY_CLUSTERS
from dark.reads import AARead

from light.database import DatabaseParameters
from light.distance import scaleLog
from light.features import Landmark
from light.landmarks import THAlphaHelix
from light.utils import stringSpans


class TestTHAlphaHelix(TestCase):
    """
    Tests for the
    light.landmark.th_alpha_helix.THAlphaHelix class.
    """

    # In the following tests, 'F' has "peak" hydrophobic, 'R' has "trough"
    # hydropathy, and 'G' has an intermediate hydropathy value.

    def check(self, sequence, expected, dbParams=None):
        """
        Check that, given a read with the sequence in C{sequence}, a
        C{THAlphaHelix} finder finds all the helices given in C{expected}.

        @param sequence: A C{str} of amino acids.
        @param expected: A C{str}, consisting of '-' and 'H' characters, where
            spans of 'H' characters indicate where helices are expected to be
            found in C{sequence}. (The non-'H' characters do not have to be
            hyphens, use whatever you like that is not an 'H'.)
        @param dbParams: A C{DatabaseParameters} instance or C{None} if default
            parameters should be used. This can be used to pass a non-default
            featureLengthBase.
        """
        self.assertEqual(len(sequence), len(expected))
        read = AARead('id', sequence)
        finder = THAlphaHelix(dbParams)
        featureLengthBase = finder._dbParams.featureLengthBase
        expectedHelices = []
        for symbol, start, end in stringSpans(expected):
            if symbol == 'H':
                length = end - start
                detail = scaleLog(length, featureLengthBase)
                expectedHelices.append(
                    Landmark('THAlphaHelix', 'THA', start, length, detail))
        self.assertEqual(expectedHelices, list(finder.find(read)))

    def testName(self):
        """
        The landmark name must be as expected.
        """
        self.assertEqual('THAlphaHelix', THAlphaHelix.NAME)

    def testSymbol(self):
        """
        The landmark symbol must be as expected.
        """
        self.assertEqual('THA', THAlphaHelix.SYMBOL)

    def testEmpty(self):
        """
        The find method must return an empty generator when given an empty
        sequence.
        """
        self.check('',
                   '')

    def testUnremarkable(self):
        """
        The find method must return an empty generator when given a sequence
        composed entirely of residues that are neither hydrophobic nor
        hydrophilic.
        """
        self.check('GGGGGGG',
                   '-------')

    def testAlternatingStartingWithHydrophobic(self):
        """
        The find method must find a helix in a sequence that alternates between
        hydrophobic and hydrophilic, starting with hydrophobic.
        """
        self.check('FRFRFRFRF',
                   'HHHHHHHHH')

    def testAlternatingStartingWithHydrophilic(self):
        """
        The find method must find a helix in a sequence that alternates between
        hydrophobic and hydrophilic, starting with hydrophilic.
        """
        self.check('RFRFRFRFR',
                   'HHHHHHHHH')

    def testNotAtStart(self):
        """
        The find method must find a helix when it is not at the start of the
        sequence.
        """
        self.check('GGGGFRRFRRFR',
                   '----HHHHHHHH')

    def testNotAtEnd(self):
        """
        The find method must find a helix when it is not at the end of the
        sequence.
        """
        self.check('FRRFRRFRGGGG',
                   'HHHHHHHH----')

    def testInMiddle(self):
        """
        The find method must find a helix in the middle of a sequence.
        """
        self.check('GGGFRRFRRFRGGGG',
                   '---HHHHHHHH----')

    def testMultipleHelices(self):
        """
        The find method must be able to find multiple helices in a sequence.
        """
        self.check('GGGFRRFRRFRGGGGGGGGFRFRFRFRFRFRGGGGGGGGFRFRFRFRGG',
                   '---HHHHHHHH--------HHHHHHHHHHHH--------HHHHHHHH--')

    def testSubsequentPeakExtendsHelixMatch(self):
        """
        If a second peak is found (with no intervening trough), it is added
        to the match.
        """
        self.check('FRFRFRFRFGGF',
                   'HHHHHHHHHHHH')

    def testSubsequentTroughExtendsHelixMatch(self):
        """
        If a second trough is found (with no intervening peak), it is added
        to the match.
        """
        self.check('FRFRFRFRGGRF',
                   'HHHHHHHHHHHH')

    def testTooShort(self):
        """
        If a potential helix is too short, it must not be found.
        """
        self.assertEqual(7, THAlphaHelix.MIN_HELIX_LENGTH)
        self.check('FRFRFR',
                   '------')

    def testInsufficientExtrema(self):
        """
        If a potential helix does not have enough extrema (total hydropathy
        peaks and troughs), it must not be found.
        """
        self.assertEqual(3, THAlphaHelix.MIN_EXTREMA_COUNT)
        self.check('FGGRGGG',
                   '-------')

    def testWaitForPeakExceeded(self):
        """
        If a potential helix does not have a peak soon enough, it must not
        be found.
        """
        self.assertEqual(6, THAlphaHelix.MAX_WAIT_FOR_PEAK)
        self.check('FGGRGGGGGF',
                   '----------')

    def testWaitForTroughExceeded(self):
        """
        If a potential helix does not have a trough soon enough, it must not
        be found.
        """
        self.assertEqual(6, THAlphaHelix.MAX_WAIT_FOR_TROUGH)
        self.check('FGGGGGRGGF',
                   '----------')

    def testPeakHydropathyCutoff(self):
        """
        The numeric value for peak hydropathy must be as expected.
        """
        self.assertEqual(3, THAlphaHelix.PEAK_HYDROPATHY)

    def testTroughHydropathyCutoff(self):
        """
        The numeric value for trough hydropathy must be as expected.
        """
        self.assertEqual(1, THAlphaHelix.TROUGH_HYDROPATHY)

    def testNonDefaultFeatureLengthBase(self):
        """
        If a non-default featureLengthBase is in the parameters given to
        the finder, it must be used.
        """
        featureLengthBase = 3.0
        dbParams = DatabaseParameters(featureLengthBase=featureLengthBase)
        read = AARead('id', 'RF' * 20)
        finder = THAlphaHelix(dbParams)
        self.assertEqual(featureLengthBase, finder._dbParams.featureLengthBase)
        # The sequence has length 40, and floor(log base 3 of 40) = 3.
        expectedSymbolDetail = scaleLog(40, featureLengthBase)
        self.assertEqual(3, expectedSymbolDetail)
        feature = list(finder.find(read))[0]
        self.assertEqual(expectedSymbolDetail, feature.symbolDetail)

    def testAlternateAAs(self):
        """
        The find method must work as expected, independent of the identity of
        the specific hydrophobic and hydrophilic residues it encounters.
        """
        peakAAs = [aa for aa in PROPERTY_CLUSTERS
                   if PROPERTY_CLUSTERS[aa]['hydropathy'] >=
                   THAlphaHelix.PEAK_HYDROPATHY]

        troughAAs = [aa for aa in PROPERTY_CLUSTERS
                     if PROPERTY_CLUSTERS[aa]['hydropathy'] <=
                     THAlphaHelix.TROUGH_HYDROPATHY]

        # Sanity check, to make sure we actually have some alternate AAs.
        self.assertEqual(7, len(peakAAs))
        self.assertEqual(7, len(troughAAs))

        for peakAA in peakAAs:
            for troughAA in troughAAs:
                self.check((peakAA + troughAA) * 4,
                           'HHHHHHHH')
