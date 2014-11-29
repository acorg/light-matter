from unittest import TestCase

from light.features import Landmark, TrigPoint


class TestLandmarks(TestCase):
    """
    Tests for the light.features.Landmark class
    """
    def testEqual(self):
        """
        Identical landmarks must compare equal.
        """
        landmark1 = Landmark('name', 'symbol', 0, 1)
        landmark2 = Landmark('name', 'symbol', 0, 1)
        self.assertEqual(landmark1, landmark2)

    def testDifferingNamesNonEqual(self):
        """
        Landmarks with different names must not compare equal.
        """
        landmark1 = Landmark('name1', 'symbol', 0, 1)
        landmark2 = Landmark('name2', 'symbol', 0, 1)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingSymbolsNonEqual(self):
        """
        Landmarks with different symbols must not compare equal.
        """
        landmark1 = Landmark('name', 'symbol1', 0, 1)
        landmark2 = Landmark('name', 'symbol2', 0, 1)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingOffsetsNonEqual(self):
        """
        Landmarks with different offsets must not compare equal.
        """
        landmark1 = Landmark('name', 'symbol', 0, 1)
        landmark2 = Landmark('name', 'symbol', 1, 1)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingLengthsNonEqual(self):
        """
        Landmarks with different lengths must not compare equal.
        """
        landmark1 = Landmark('name', 'symbol1', 0, 1)
        landmark2 = Landmark('name', 'symbol2', 0, 2)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingRepeatCountssNonEqual(self):
        """
        Landmarks with different repeat counts must not compare equal.
        """
        landmark1 = Landmark('name', 'symbol1', 0, 1, 0)
        landmark2 = Landmark('name', 'symbol2', 0, 1, 1)
        self.assertNotEqual(landmark1, landmark2)


class TestTrigPoints(TestCase):
    """
    Tests for the light.features.TrigPoint class
    """
    def testEqual(self):
        """
        Identical landmarks must compare equal.
        """
        landmark1 = TrigPoint('name', 'symbol', 0)
        landmark2 = TrigPoint('name', 'symbol', 0)
        self.assertEqual(landmark1, landmark2)

    def testDifferingNamesNonEqual(self):
        """
        Trig points with different names must not compare equal.
        """
        landmark1 = TrigPoint('name1', 'symbol', 0)
        landmark2 = TrigPoint('name2', 'symbol', 0)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingSymbolsNonEqual(self):
        """
        Trig points with different symbols must not compare equal.
        """
        landmark1 = TrigPoint('name', 'symbol1', 0)
        landmark2 = TrigPoint('name', 'symbol2', 0)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingOffsetsNonEqual(self):
        """
        Trig points with different offsets must not compare equal.
        """
        landmark1 = TrigPoint('name', 'symbol', 0)
        landmark2 = TrigPoint('name', 'symbol', 1)
        self.assertNotEqual(landmark1, landmark2)

    def testDifferingLengthsNonEqual(self):
        """
        Trig points with different lengths must not compare equal.
        """
        landmark1 = TrigPoint('name', 'symbol1', 0)
        landmark2 = TrigPoint('name', 'symbol2', 0)
        self.assertNotEqual(landmark1, landmark2)
