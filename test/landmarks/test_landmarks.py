from unittest import TestCase

from light.landmarks import (
    findLandmark, ALL_LANDMARK_FINDER_CLASSES, DEFAULT_LANDMARK_FINDER_CLASSES)
from light.landmarks.alpha_helix import AlphaHelix
from light.landmarks.alpha_helix_3_10 import AlphaHelix_3_10
from light.landmarks.alpha_helix_pi import AlphaHelix_pi
from light.landmarks.beta_strand import BetaStrand
from light.landmarks.turn import Turn


class TestFindLandmark(TestCase):
    """
    Tests for the light.landmarks.findLandmark function.
    """

    def testFindFails(self):
        """
        The find function should return C{None} if asked to find a landmark
        class that doesn't exist.
        """
        self.assertIs(None, findLandmark('silly'))

    def testFindAllClasses(self):
        """
        The find function should be able to find all landmark classes by name.
        """
        for klass in ALL_LANDMARK_FINDER_CLASSES:
            self.assertIs(klass, findLandmark(klass.NAME))


class TestAllLandmarkClasses(TestCase):
    """
    Trivial tests for the ALL_LANDMARK_FINDER_CLASSES set.
    """

    def testAllClasses(self):
        """
        The ALL_LANDMARK_FINDER_CLASSES set must be as expected.
        """
        self.assertEqual(
            {AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand, Turn},
            ALL_LANDMARK_FINDER_CLASSES)


class TestDefaultLandmarkClasses(TestCase):
    """
    Trivial tests for the DEFAULT_LANDMARK_FINDER_CLASSES set.
    """

    def testDefaultClasses(self):
        """
        The DEFAULT_LANDMARK_FINDER_CLASSES must be as expected.
        """
        self.assertEqual(
            {AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand, Turn},
            DEFAULT_LANDMARK_FINDER_CLASSES)

    def testDefaultClassesAreInAllClasses(self):
        """
        The DEFAULT_LANDMARK_FINDER_CLASSES must all appear in
        ALL_LANDMARK_FINDER_CLASSES.
        """
        for klass in DEFAULT_LANDMARK_FINDER_CLASSES:
            self.assertIn(klass, ALL_LANDMARK_FINDER_CLASSES)
