from unittest import TestCase

from dark.reads import AARead

from light.distance import scale
from light.features import Landmark
from light.parameters import Parameters
from light.landmarks.gor4_coil import GOR4Coil


class TestGOR4Coil(TestCase):
    """
    Tests for the landmark.gor4_coil.GOR4Coil class.
    """
    def testClassAttributes(self):
        """
        The GOR4Coil class attributes must be as expected.
        """
        self.assertEqual('GOR4Coil', GOR4Coil.NAME)
        self.assertEqual('GC', GOR4Coil.SYMBOL)

    def testCoilInMiddleAndAtEnd(self):
        """
        The GOR4Coil landmark finder must find a coil in the middle and at the
        end of a sequence.
        """
        read = AARead('id', 'VICVIC')
        landmark = GOR4Coil()
        result = list(landmark.find(read))
        len1 = scale(1, Parameters.DEFAULT_DISTANCE_BASE)
        # The GOR IV secondary structure prediction is 'EECEEC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 2, 1, len1),
                          Landmark('GOR4Coil', 'GC', 5, 1, len1)],
                         result)

    def testAllCoil(self):
        """
        The GOR4Coil landmark finder must find a coil that spans the whole
        sequence.
        """
        read = AARead('id', 'EA')
        landmark = GOR4Coil()
        result = list(landmark.find(read))
        len2 = scale(2, Parameters.DEFAULT_DISTANCE_BASE)
        # The GOR IV secondary structure prediction is 'CC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, len2, len2)],
                         result)

    def testApoamicyaninFiveCoils(self):
        """
        The GOR4Coil landmark finder must find the five expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database.
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        landmark = GOR4Coil()
        result = list(landmark.find(read))
        len1 = scale(1, Parameters.DEFAULT_DISTANCE_BASE)
        len2 = scale(2, Parameters.DEFAULT_DISTANCE_BASE)
        len4 = scale(4, Parameters.DEFAULT_DISTANCE_BASE)
        len10 = scale(10, Parameters.DEFAULT_DISTANCE_BASE)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 10, len10),
                          Landmark('GOR4Coil', 'GC', 17, 2, len2),
                          Landmark('GOR4Coil', 'GC', 30, 4, len4),
                          Landmark('GOR4Coil', 'GC', 39, 2, len2),
                          Landmark('GOR4Coil', 'GC', 49, 1, len1)],
                         result)

    def testApoamicyaninFiveCoilsWithNonDefaultDistanceBase(self):
        """
        The GOR4Coil landmark finder must find the five expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database. It must return the right length of the landmark after a
        non-default distanceBase has been applied.
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        distanceBase = 1.5
        landmark = GOR4Coil(distanceBase=distanceBase)
        result = list(landmark.find(read))
        len1 = scale(1, distanceBase)
        len2 = scale(2, distanceBase)
        len4 = scale(4, distanceBase)
        len10 = scale(10, distanceBase)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'.
        self.maxDiff = None
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 10, len10),
                          Landmark('GOR4Coil', 'GC', 17, 2, len2),
                          Landmark('GOR4Coil', 'GC', 30, 4, len4),
                          Landmark('GOR4Coil', 'GC', 39, 2, len2),
                          Landmark('GOR4Coil', 'GC', 49, 1, len1)],
                         result)

    def testLengthMustBeStoredCorrectly(self):
        """
        The length of a GOR4BetaStrand must be stored correctly (not scaled).
        """
        read = AARead('id', 'DKATIPSESP')
        landmark = GOR4Coil()
        result = list(landmark.find(read))
        len6 = scale(6, Parameters.DEFAULT_DISTANCE_BASE)
        len1 = scale(1, Parameters.DEFAULT_DISTANCE_BASE)
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 6, len6),
                          Landmark('GOR4Coil', 'GC', 9, 1, len1)],
                         result)
