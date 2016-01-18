from unittest import TestCase

from dark.reads import AARead

from light.landmarks.random import RandomLandmark


class TestRandomLandmark(TestCase):
    """
    Tests for the Landmark.RandomLandmark class.
    """
    def testCorrectNumberOfLandmarksDefaultDensity(self):
        """
        The correct number of landmarks must be returned, as specified by the
        default density parameter.
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
        landmark = RandomLandmark()
        result = list(landmark.find(read, density=0.5))
        self.assertEqual(10, len(result))

    def testAttributes(self):
        """
        The returned landmarks must have the correct attributes.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        landmark = RandomLandmark()
        result = list(landmark.find(read, density=0.05))[0]
        self.assertEqual('RandomLandmark', result.name)
        self.assertEqual('RL', result.symbol)
        self.assertEqual(1, result.length)
        self.assertTrue(0 <= result.offset < len(read))
