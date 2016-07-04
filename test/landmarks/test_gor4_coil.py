from unittest import TestCase

from dark.reads import AARead

from light.distance import scaleLog
from light.features import Landmark
from light.landmarks.gor4_coil import GOR4Coil
from light.parameters import DatabaseParameters


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
        scaled2 = scaleLog(2, DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE)
        scaled0 = scaleLog(0, DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE)
        # The GOR IV secondary structure prediction is 'EECEEC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 2, scaled2),
                          Landmark('GOR4Coil', 'GC', 5, 1, scaled0)],
                         result)

    def testAllCoil(self):
        """
        The GOR4Coil landmark finder must find a coil that spans the whole
        sequence.
        """
        read = AARead('id', 'EA')
        landmark = GOR4Coil()
        result = list(landmark.find(read))
        scaled2 = scaleLog(2, DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE)
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
        landmark = GOR4Coil()
        result = list(landmark.find(read))
        scaled2 = scaleLog(2, DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE)
        scaled3 = scaleLog(3, DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE)
        scaled5 = scaleLog(5, DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE)
        scaled10 = scaleLog(10, DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 10, scaled10),
                          Landmark('GOR4Coil', 'GC', 17, 2, scaled2),
                          Landmark('GOR4Coil', 'GC', 28, 6, scaled5),
                          Landmark('GOR4Coil', 'GC', 39, 3, scaled3),
                          Landmark('GOR4Coil', 'GC', 45, 5, scaled5)],
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
        dbParams = DatabaseParameters(featureLengthBase=featureLengthBase)
        landmark = GOR4Coil(dbParams)
        result = list(landmark.find(read))
        scaled2 = scaleLog(2, featureLengthBase)
        scaled3 = scaleLog(3, featureLengthBase)
        scaled4 = scaleLog(4, featureLengthBase)
        scaled6 = scaleLog(6, featureLengthBase)
        scaled10 = scaleLog(10, featureLengthBase)
        # The GOR IV secondary structure prediction is
        # 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'.
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 10, scaled10),
                          Landmark('GOR4Coil', 'GC', 17, 2, scaled2),
                          Landmark('GOR4Coil', 'GC', 28, 6, scaled6),
                          Landmark('GOR4Coil', 'GC', 39, 3, scaled3),
                          Landmark('GOR4Coil', 'GC', 45, 5, scaled4)],
                         result)

    def testLengthMustBeStoredCorrectly(self):
        """
        The length of a GOR4BetaStrand must be stored correctly (not scaled).
        """
        read = AARead('id', 'DKATIPSESP')
        landmark = GOR4Coil()
        result = list(landmark.find(read))
        scaled9 = scaleLog(9, DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE)
        self.assertEqual([Landmark('GOR4Coil', 'GC', 0, 10, scaled9)],
                         result)
