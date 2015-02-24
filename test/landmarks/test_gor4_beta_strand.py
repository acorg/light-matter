from unittest import TestCase

from dark.reads import AARead

from light.distance import scale
from light.features import Landmark
from light.database import Database
from light.landmarks.gor4_beta_strand import GOR4BetaStrand


class TestGOR4BetaStrand(TestCase):
    """
    Tests for the landmark.gor4_beta_strand.GOR4BetaStrand class.
    """
    def testClassAttributes(self):
        """
        The GOR4BetaStrand class attributes must be as expected.
        """
        self.assertEqual('GOR4BetaStrand', GOR4BetaStrand.NAME)
        self.assertEqual('GB', GOR4BetaStrand.SYMBOL)

    def testNoBetaStrands(self):
        """
        The GOR4BetaStrand landmark finder must not find any beta strands when
        none are present.
        """
        read = AARead('id', 'EA')
        landmark = GOR4BetaStrand()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testBetaStrandAtStartAndInMiddle(self):
        """
        The GOR4BetaStrand landmark finder must find a beta strand that
        occurs at the start of a sequence, as well as one that appears
        in the middle.
        """
        read = AARead('id', 'VICVICV')
        landmark = GOR4BetaStrand()
        result = list(landmark.find(read))
        len2 = scale(2, Database.DEFAULT_DISTANCE_BASE)
        len3 = scale(3, Database.DEFAULT_DISTANCE_BASE)
        # The GOR IV secondary structure prediction is 'EECEEEC'.
        self.assertEqual([Landmark('GOR4BetaStrand', 'GB', 0, len2, len2),
                          Landmark('GOR4BetaStrand', 'GB', 3, len3, len3)],
                         result)

    def testBetaStrandAtEnd(self):
        """
        The GOR4BetaStrand landmark finder must find a beta strand that
        occurs at the end of a sequence.
        """
        read = AARead('id', 'LHV')
        landmark = GOR4BetaStrand()
        result = list(landmark.find(read))
        len2 = scale(2, Database.DEFAULT_DISTANCE_BASE)
        # The GOR IV secondary structure prediction is 'CEE'.
        self.assertEqual([Landmark('GOR4BetaStrand', 'GB', 1, len2, len2)],
                         result)

    def testApoamicyaninTwoBetaStrands(self):
        """
        The GOR4BetaStrand landmark finder must find the two expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database.
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        landmark = GOR4BetaStrand()
        result = list(landmark.find(read))
        len5 = scale(5, Database.DEFAULT_DISTANCE_BASE)
        len8 = scale(8, Database.DEFAULT_DISTANCE_BASE)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'
        self.assertEqual([Landmark('GOR4BetaStrand', 'GB', 34, len5, len5),
                          Landmark('GOR4BetaStrand', 'GB', 41, len8, len8)],
                         result)

    def testApoamicyaninTwoBetaStrandsWithNonDefaultDistanceBase(self):
        """
        The GOR4BetaStrand landmark finder must find the two expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database. It must return the right length of the landmark, after a
        non-default distanceBase was applied.
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        distanceBase = 1.5
        landmark = GOR4BetaStrand(distanceBase=distanceBase)
        result = list(landmark.find(read))
        len5 = scale(5, distanceBase)
        len8 = scale(8, distanceBase)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'
        self.assertEqual([Landmark('GOR4BetaStrand', 'GB', 34, len5, len5),
                          Landmark('GOR4BetaStrand', 'GB', 41, len8, len8)],
                         result)
