from unittest import TestCase

from light.landmarks import find, ALL_LANDMARK_FINDER_CLASSES
from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.alpha_helix_3_10 import AlphaHelix_3_10
from light.landmarks.alpha_helix_pi import AlphaHelix_pi


class TestFind(TestCase):
    """
    Tests for the light.landmarks.find function.
    """

    def testFindFails(self):
        """
        The find function should return C{None} if asked to find a landmark
        class that doesn't exist.
        """
        self.assertIs(None, find('silly'))

    def testFindAlphaHelix(self):
        """
        The find function should be able to find the AlphaHelix class by name.
        """
        self.assertIs(AlphaHelix, find('AlphaHelix'))


class TestAllClasses(TestCase):
    """
    Trivial tests for the ALL_LANDMARK_FINDER_CLASSES set.
    """

    def testAllClassesContainsAlphaHelix(self):
        """
        The ALL_LANDMARK_FINDER_CLASSES set must contain the AlphaHelix class.
        """
        self.assertIn(AlphaHelix, ALL_LANDMARK_FINDER_CLASSES)

    def testAllClassesContainsAlphaHelix_3_10(self):
        """
        The ALL_LANDMARK_FINDER_CLASSES set must contain the AlphaHelix_3_10
        class.
        """
        self.assertIn(AlphaHelix_3_10, ALL_LANDMARK_FINDER_CLASSES)

    def testAllClassesContainsAlphaHelix_pi(self):
        """
        The ALL_LANDMARK_FINDER_CLASSES set must contain the AlphaHelix_pi
        class.
        """
        self.assertIn(AlphaHelix_pi, ALL_LANDMARK_FINDER_CLASSES)
