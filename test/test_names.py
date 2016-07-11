from unittest import TestCase

from light.landmarks import ALL_LANDMARK_CLASSES_INCLUDING_DEV
from light.trig import ALL_TRIG_CLASSES_INCLUDING_DEV


class TestNames(TestCase):
    """
    Tests to ensure that all landmark and trig point finders have distinct
    names.
    """
    def testRightNumberOfNames(self):
        """
        Must find the right number of names.
        """
        trigNames = [name.NAME for name in ALL_TRIG_CLASSES_INCLUDING_DEV]
        lmNames = [name.NAME for name in ALL_LANDMARK_CLASSES_INCLUDING_DEV]
        names = trigNames + lmNames
        self.assertEqual(29, len(names))

    def testAllLandmarkNamesAreDistinct(self):
        """
        All landmark finder names must be distinct.
        """
        names = [name.NAME for name in ALL_LANDMARK_CLASSES_INCLUDING_DEV]
        self.assertEqual(len(names), len(set(names)))

    def testAllTrigNamesAreDistinct(self):
        """
        All trig point finder names must be distinct.
        """
        names = [name.NAME for name in ALL_TRIG_CLASSES_INCLUDING_DEV]
        self.assertEqual(len(names), len(set(names)))

    def testLandmarkAndTrigNamesAreDistinct(self):
        """
        All landmark and trig point finder names must be distinct.
        """
        trigNames = [name.NAME for name in ALL_TRIG_CLASSES_INCLUDING_DEV]
        lmNames = [name.NAME for name in ALL_LANDMARK_CLASSES_INCLUDING_DEV]
        names = trigNames + lmNames
        self.assertEqual(len(names), len(set(names)))
