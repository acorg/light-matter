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
        scaled1 = scale(1, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        # The GOR IV secondary structure prediction is 'EECEEC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 2, 1, scaled1),
                          Landmark('GOR4Coil', 'GC', 5, 1, scaled1)],
                         result)

    def testAllCoil(self):
        """
        The GOR4Coil landmark finder must find a coil that spans the whole
        sequence.
        """
        read = AARead('id', 'EA')
        landmark = GOR4Coil()
        result = list(landmark.find(read))
        scaled2 = scale(2, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        # The GOR IV secondary structure prediction is 'CC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 2, scaled2)],
                         result)

    def testApoamicyaninFiveCoils(self):
        """
        The GOR4Coil landmark finder must find the five expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database.
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        landmark = GOR4Coil(Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        result = list(landmark.find(read))
        scaled1 = scale(1, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        scaled2 = scale(2, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        scaled4 = scale(4, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        scaled10 = scale(10, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 10, scaled10),
                          Landmark('GOR4Coil', 'GC', 17, 2, scaled2),
                          Landmark('GOR4Coil', 'GC', 30, 4, scaled4),
                          Landmark('GOR4Coil', 'GC', 39, 2, scaled2),
                          Landmark('GOR4Coil', 'GC', 49, 1, scaled1)],
                         result)

    def testApoamicyaninFiveCoilsWithNonDefaultFeatureLengthBase(self):
        """
        The GOR4Coil landmark finder must find the five expected landmarks
        in a fragment of the APOAMICYANIN sequence from the GOR IV reference
        database. It must have the right scaled length of the landmark after a
        non-default featureLengthBase has been applied.
        """
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        featureLengthBase = 1.5
        landmark = GOR4Coil(featureLengthBase)
        result = list(landmark.find(read))
        scaled1 = scale(1, featureLengthBase)
        scaled2 = scale(2, featureLengthBase)
        scaled4 = scale(4, featureLengthBase)
        scaled10 = scale(10, featureLengthBase)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 10, scaled10),
                          Landmark('GOR4Coil', 'GC', 17, 2, scaled2),
                          Landmark('GOR4Coil', 'GC', 30, 4, scaled4),
                          Landmark('GOR4Coil', 'GC', 39, 2, scaled2),
                          Landmark('GOR4Coil', 'GC', 49, 1, scaled1)],
                         result)

    def testLengthMustBeStoredCorrectly(self):
        """
        The length of a GOR4BetaStrand must be stored correctly (not scaled).
        """
        read = AARead('id', 'DKATIPSESP')
        landmark = GOR4Coil()
        result = list(landmark.find(read))
        scaled6 = scale(6, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        scaled1 = scale(1, Parameters.DEFAULT_FEATURE_LENGTH_BASE)
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 6, scaled6),
                          Landmark('GOR4Coil', 'GC', 9, 1, scaled1)],
                         result)
