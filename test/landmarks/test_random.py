from unittest import TestCase

from dark.reads import AARead

from light.landmarks.random import RandomLandmark
from light.parameters import DatabaseParameters


class TestRandomLandmark(TestCase):
    """
    Tests for the Landmark.RandomLandmark class.
    """
    def testCorrectNumberOfLandmarksDefaultDensity(self):
        """
        The correct number of landmarks must be returned, based on the
        default density parameter value (0.1).
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        landmark = RandomLandmark()
        result = list(landmark.find(read))
        self.assertEqual(2, len(result))

    def testCorrectNumberOfLandmarksNonDefaultDensity(self):
        """
        The correct number of landmarks must be returned, as specified by the
        density parameter.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        dbParams = DatabaseParameters(randomLandmarkDensity=0.5)
        landmark = RandomLandmark(dbParams)
        result = list(landmark.find(read))
        self.assertEqual(10, len(result))

    def testAttributes(self):
        """
        The returned landmarks must have the correct attributes.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        # Use an extremely high density, to make sure the finder never
        # fails to produce a random landmark.
        dbParams = DatabaseParameters(randomLandmarkDensity=0.9)
        landmark = RandomLandmark(dbParams)
        result = list(landmark.find(read))[0]
        self.assertEqual('RandomLandmark', result.name)
        self.assertEqual('RL', result.symbol)
        self.assertEqual(1, result.length)
        self.assertTrue(0 <= result.offset < len(read))
