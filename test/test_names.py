from unittest import TestCase

from light.landmarks import ALL_LANDMARK_CLASSES
from light.trig import ALL_TRIG_CLASSES


class TestNames(TestCase):
    """
    Tests to ensure that all landmark and trig point finders have different
    names.
    """
    def testRightNumberOfNames(self):
        """
        Must find the right number of names.
        """
        trigNames = [name.NAME for name in ALL_TRIG_CLASSES]
        lmNames = [name.NAME for name in ALL_LANDMARK_CLASSES]
        names = trigNames + lmNames
        self.assertEqual(16, len(names))

    def testAllLandmarkNamesAreDifferent(self):
        """
        All landmark finder names must be different.
        """
        names = [name.NAME for name in ALL_LANDMARK_CLASSES]
        self.assertEqual(len(names), len(set(names)))

    def testAllTrigNamesAreDifferent(self):
        """
        All trig point finder names must be different.
        """
        names = [name.NAME for name in ALL_TRIG_CLASSES]
        self.assertEqual(len(names), len(set(names)))

    def testLandmarkAndTrigNamesAreDifferent(self):
        """
        All landmark and trig point finder names must be different.
        """
        trigNames = [name.NAME for name in ALL_TRIG_CLASSES]
        lmNames = [name.NAME for name in ALL_LANDMARK_CLASSES]
        names = trigNames + lmNames
        self.assertEqual(len(names), len(set(names)))
