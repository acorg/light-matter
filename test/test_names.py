from unittest import TestCase

from light.landmarks import ALL_LANDMARK_FINDER_CLASSES
from light.trig import ALL_TRIG_FINDER_CLASSES


class TestNames(TestCase):
    """
    Tests to ensure that all landmarks and trig points have different names.
    """
    def testRightNumberOfNames(self):
        """
        Must find the right number of names.
        """
        trigNames = [name.NAME for name in ALL_TRIG_FINDER_CLASSES]
        lmNames = [name.NAME for name in ALL_LANDMARK_FINDER_CLASSES]
        names = trigNames + lmNames
        self.assertEqual(10, len(names))

    def testAllLandmarkNamesAreDifferent(self):
        """
        All landmark names must be different.
        """
        names = [name.NAME for name in ALL_LANDMARK_FINDER_CLASSES]
        self.assertEqual(len(names), len(set(names)))

    def testAllTrigNamesAreDifferent(self):
        """
        All trig finder names must be different.
        """
        names = [name.NAME for name in ALL_TRIG_FINDER_CLASSES]
        self.assertEqual(len(names), len(set(names)))

    def testLandmarkAndTrigNamesAreDifferent(self):
        """
        All landmark and trig finder names must be different.
        """
        trigNames = [name.NAME for name in ALL_TRIG_FINDER_CLASSES]
        lmNames = [name.NAME for name in ALL_LANDMARK_FINDER_CLASSES]
        names = trigNames + lmNames
        self.assertEqual(len(names), len(set(names)))
