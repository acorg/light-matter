from unittest import TestCase

from light.landmarks import (
    findLandmark, findLandmarks, ALL_LANDMARK_CLASSES,
    DEFAULT_LANDMARK_CLASSES, AlphaHelix, AlphaHelix_3_10,
    AlphaHelix_pi, AminoAcids, BetaStrand, BetaTurn, GOR4AlphaHelix,
    GOR4BetaStrand, GOR4Coil, Prosite)


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
        for klass in ALL_LANDMARK_CLASSES:
            self.assertIs(klass, findLandmark(klass.NAME))


class TestFindLandmarks(TestCase):
    """
    Tests for the light.landmarks.findLandmarks function.
    """

    def testFindFails(self):
        """
        The find function should raise a ValueError if asked to find unknown
        classes.
        """
        error = '^Unknown landmark finders: x, y\.$'
        self.assertRaisesRegexp(ValueError, error, findLandmarks, ['x', 'y'])

    def testFindNone(self):
        """
        The find function should return an empty list if passed None.
        """
        result = findLandmarks(None)
        self.assertEqual([], result)

    def testFindKnownClasses(self):
        """
        The find function should return known classes correctly.
        """
        result = findLandmarks(['AlphaHelix', 'BetaStrand'])
        self.assertEqual(2, len(result))
        self.assertIs(AlphaHelix, result[0])
        self.assertIs(BetaStrand, result[1])


class TestAllLandmarkClasses(TestCase):
    """
    Trivial tests for the ALL_LANDMARK_CLASSES set.
    """

    def testAllClasses(self):
        """
        The ALL_LANDMARK_CLASSES set must be as expected.
        """
        self.assertEqual(
            {AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, AminoAcids,
             BetaStrand, BetaTurn, GOR4AlphaHelix, GOR4BetaStrand, GOR4Coil,
             Prosite},
            ALL_LANDMARK_CLASSES)


class TestDefaultLandmarkClasses(TestCase):
    """
    Trivial tests for the DEFAULT_LANDMARK_CLASSES set.
    """

    def testDefaultClasses(self):
        """
        The DEFAULT_LANDMARK_CLASSES must be as expected.
        """
        self.assertEqual(
            {AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, AminoAcids,
             BetaStrand, BetaTurn, Prosite},
            DEFAULT_LANDMARK_CLASSES)

    def testDefaultClassesAreInAllClasses(self):
        """
        The DEFAULT_LANDMARK_CLASSES must all appear in
        ALL_LANDMARK_CLASSES.
        """
        for klass in DEFAULT_LANDMARK_CLASSES:
            self.assertIn(klass, ALL_LANDMARK_CLASSES)
