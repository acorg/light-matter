from unittest import TestCase

from dark.reads import AARead

from light.distance import scale
from light.features import Landmark
from light.parameters import Parameters
from light.database import Database
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
        len7 = scale(7, Parameters.DEFAULT_DISTANCE_BASE)
        len11 = scale(11, Parameters.DEFAULT_DISTANCE_BASE)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'
        self.assertEqual([Landmark('GOR4AlphaHelix', 'GA', 10, len7, len7),
                          Landmark('GOR4AlphaHelix', 'GA', 19, len11, len11)],
                         result)

    def testApoamicyaninTwoAlphaHelixsWithNonDefaultDistanceBase(self):
        """
        The GOR4AlphaHelix landmark finder must find the two expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database. It must return the right length of the landmark, after a
        non-default distanceBase has been applied
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        distanceBase = 1.5
        landmark = GOR4AlphaHelix(distanceBase=distanceBase)
        result = list(landmark.find(read))
        len7 = scale(7, distanceBase)
        len11 = scale(11, distanceBase)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'
        self.assertEqual([Landmark('GOR4AlphaHelix', 'GA', 10, len7, len7),
                          Landmark('GOR4AlphaHelix', 'GA', 19, len11, len11)],
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
        alphaHelixDb = Database(Parameters([GOR4AlphaHelix], []))
        betaStrandDb = Database(Parameters([GOR4BetaStrand], []))
        alphaHelixScanned = alphaHelixDb.scan(self.READ)
        betaStrandScanned = betaStrandDb.scan(self.READ)
        alphaHelixIndices = alphaHelixScanned.coveredIndices()
        betaStrandIndices = betaStrandScanned.coveredIndices()
        self.assertEqual(0, len(alphaHelixIndices & betaStrandIndices))

    def testNoOverlapDistanceBaseOne(self):
        """
        There cannot be any index overlap between landmarks found by the
        GOR4 alpha helix and beta strand finders using a distance base of 1.0
        (which should do no scaling).
        """
        alphaHelixDb = Database(
            Parameters([GOR4AlphaHelix], [], distanceBase=1.0))
        betaStrandDb = Database(
            Parameters([GOR4BetaStrand], [], distanceBase=1.0))
        alphaHelixScanned = alphaHelixDb.scan(self.READ)
        betaStrandScanned = betaStrandDb.scan(self.READ)
        alphaHelixIndices = alphaHelixScanned.coveredIndices()
        betaStrandIndices = betaStrandScanned.coveredIndices()
        self.assertEqual(0, len(alphaHelixIndices & betaStrandIndices))

    def testNoOverlapDistanceBaseOnePointFive(self):
        """
        There cannot be any index overlap between landmarks found by the
        GOR4 alpha helix and beta strand finders using a distance base of 1.5.
        """
        alphaHelixDb = Database(
            Parameters([GOR4AlphaHelix], [], distanceBase=1.5))
        betaStrandDb = Database(
            Parameters([GOR4BetaStrand], [], distanceBase=1.5))
        alphaHelixScanned = alphaHelixDb.scan(self.READ)
        betaStrandScanned = betaStrandDb.scan(self.READ)
        alphaHelixIndices = alphaHelixScanned.coveredIndices()
        betaStrandIndices = betaStrandScanned.coveredIndices()
        self.assertEqual(0, len(alphaHelixIndices & betaStrandIndices))
