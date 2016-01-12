from unittest import TestCase

from dark.reads import AARead

from light.distance import scale
from light.features import Landmark
from light.parameters import Parameters
from light.backend import Backend
from light.landmarks.gor4_alpha_helix import GOR4AlphaHelix
from light.landmarks.gor4_beta_strand import GOR4BetaStrand


class TestGOR4AlphaHelix(TestCase):
    """
    Tests for the landmark.gor4_alpha_helix.GOR4AlphaHelix class.
    """
    def testClassAttributes(self):
        """
        The GOR4AlphaHelix class attributes must be as expected.
        """
        self.assertEqual('GOR4AlphaHelix', GOR4AlphaHelix.NAME)
        self.assertEqual('GA', GOR4AlphaHelix.SYMBOL)

    def testNoAlphaHelices(self):
        """
        The GOR4AlphaHelix landmark finder must not find any alpha helices when
        none are present.
        """
        read = AARead('id', 'EA')
        landmark = GOR4AlphaHelix()
        result = list(landmark.find(read))
        self.assertEqual([], result)

    def testApoamicyaninTwoAlphaHelixs(self):
        """
        The GOR4AlphaHelix landmark finder must find the two expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database.
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        landmark = GOR4AlphaHelix()
        result = list(landmark.find(read))
        scaled7 = scale(7, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        scaled11 = scale(11, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'
        self.assertEqual([Landmark('GOR4AlphaHelix', 'GA', 10, 7, scaled7),
                          Landmark('GOR4AlphaHelix', 'GA', 19, 11, scaled11)],
                         result)

    def testApoamicyaninTwoAlphaHelixsWithNonDefaultFeatureLengthBase(self):
        """
        The GOR4AlphaHelix landmark finder must find the two expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database. It must have the right length of the landmark, after a
        non-default featureLengthBase has been applied
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        featureLengthBase = 1.5
        landmark = GOR4AlphaHelix(featureLengthBase)
        result = list(landmark.find(read))
        scaled7 = scale(7, featureLengthBase)
        scaled11 = scale(11, featureLengthBase)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'
        self.assertEqual([Landmark('GOR4AlphaHelix', 'GA', 10, 7, scaled7),
                          Landmark('GOR4AlphaHelix', 'GA', 19, 11, scaled11)],
                         result)

    def testStoreLengthCorrectly(self):
        """
        The length of a GOR4BetaStrand must be stored correctly (not scaled).
        """
        seq = 'DKATIPSESPFAAAEVAAIVFAAAEVAAIVVFAAAEVAAIVVDIAKMKYFAAAEVAAIVVDI'
        read = AARead('id', seq)
        featureLengthBase = 1.5
        landmark = GOR4AlphaHelix(featureLengthBase)
        result = list(landmark.find(read))
        scaled47 = scale(47, featureLengthBase)
        self.assertEqual([Landmark('GOR4AlphaHelix', 'GA', 10, 47, scaled47)],
                         result)


class TestGOR4BetaStrandOverlap(TestCase):
    """
    Test that the GOR4AlphaHelix finder does not produce landmarks that
    overlap with those of the GOR4BetaStrand finder.

    This is for debugging https://github.com/acorg/light-matter/issues/220
    """
    READ = AARead('gi|188036137|pdb|2VQ0_A|Chain A, Capsid Structure Of '
                  'Sesbania Mosaic Virus Coat Protein Deletion Mutant '
                  'Rcp(Delta 48 To 59)',
                  'MAKRLSKQQLAKAIANTLETPPQPKAGRRRNRRRQRSAVQQLQPTQAVRIRNP'
                  'AVSSSRGGITVLTHSELSAEIGVTDSIVVSSELVMPYTVGTWLRGVAANWSKY'
                  'SWLSVRYTYIPSCPSSTAGSIHMGFQYDMADTVPVSVNQLSNLRGYVSGQVWS'
                  'GSAGLCFINGTRCSDTSTAISTTLDVSKLGKKWYPYKTSADYATAVGVDVNIA'
                  'TPLVPARLVIALLDGSSSTAVAAGRIYCTYTIQMIEPTASALNN')

    def testNoOverlapDefaultDistanceBase(self):
        """
        There cannot be any index overlap between landmarks found by the
        GOR4 alpha helix and beta strand finders using the default distance
        base (currently 1.1).
        """
        alphaHelixBe = Backend()
        alphaHelixBe.configure(Parameters([GOR4AlphaHelix], []))
        betaStrandBe = Backend()
        betaStrandBe.configure(Parameters([GOR4BetaStrand], []))
        alphaHelixScanned = alphaHelixBe.scan(self.READ)
        betaStrandScanned = betaStrandBe.scan(self.READ)
        alphaHelixIndices = alphaHelixScanned.coveredIndices()
        betaStrandIndices = betaStrandScanned.coveredIndices()
        self.assertEqual(0, len(alphaHelixIndices & betaStrandIndices))

    def testNoOverlapDistanceBaseOne(self):
        """
        There cannot be any index overlap between landmarks found by the
        GOR4 alpha helix and beta strand finders using a distance base of 1.0
        (which should do no scaling).
        """
        alphaHelixBe = Backend()
        alphaHelixBe.configure(
            Parameters([GOR4AlphaHelix], [], distanceBase=1.0))
        betaStrandBe = Backend()
        betaStrandBe.configure(
            Parameters([GOR4BetaStrand], [], distanceBase=1.0))
        alphaHelixScanned = alphaHelixBe.scan(self.READ)
        betaStrandScanned = betaStrandBe.scan(self.READ)
        alphaHelixIndices = alphaHelixScanned.coveredIndices()
        betaStrandIndices = betaStrandScanned.coveredIndices()
        self.assertEqual(0, len(alphaHelixIndices & betaStrandIndices))

    def testNoOverlapDistanceBaseOnePointFive(self):
        """
        There cannot be any index overlap between landmarks found by the
        GOR4 alpha helix and beta strand finders using a distance base of 1.5.
        """
        alphaHelixBe = Backend()
        alphaHelixBe.configure(
            Parameters([GOR4AlphaHelix], [], distanceBase=1.5))
        betaStrandBe = Backend()
        betaStrandBe.configure(
            Parameters([GOR4BetaStrand], [], distanceBase=1.5))
        alphaHelixScanned = alphaHelixBe.scan(self.READ)
        betaStrandScanned = betaStrandBe.scan(self.READ)
        alphaHelixIndices = alphaHelixScanned.coveredIndices()
        betaStrandIndices = betaStrandScanned.coveredIndices()
        self.assertEqual(0, len(alphaHelixIndices & betaStrandIndices))
