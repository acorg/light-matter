from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.landmarks.gor4_alpha_helix import GOR4AlphaHelix


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
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'
        self.assertEqual([Landmark('GOR4AlphaHelix', 'GA', 10, 7, 7),
                          Landmark('GOR4AlphaHelix', 'GA', 19, 11, 11)],
                         result)

    def testApoamicyaninTwoAlphaHelixsWithBucketFactor(self):
        """
        The GOR4AlphaHelix landmark finder must find the two expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database. It must return the right length of the landmark, after a
        distanceScale has been applied
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        landmark = GOR4AlphaHelix(distanceScale=1.5)
        result = list(landmark.find(read))
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'
        self.assertEqual([Landmark('GOR4AlphaHelix', 'GA', 10, 4, 4),
                          Landmark('GOR4AlphaHelix', 'GA', 19, 7, 7)],
                         result)
