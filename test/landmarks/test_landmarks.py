from unittest import TestCase

from light.landmarks import (
    findLandmark, findLandmarks, landmarkNameFromHashkey,
    ALL_LANDMARK_CLASSES, DEFAULT_LANDMARK_CLASSES, AlphaHelix,
    AlphaHelix_3_10, AlphaHelix_pi, AminoAcids, BetaStrand, BetaTurn,
    GOR4AlphaHelix, GOR4BetaStrand, GOR4Coil, Prosite)


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


class TestLandmarkNameFromHashkey(TestCase):
    """
    Tests for the light.landmarks.landmarkNameFromHashkey function.
    """

    def testFail(self):
        """
        The landmarkNameFromHashkey function should return None if asked to
        find a hashkey that no class created.
        """
        self.assertIs(None, landmarkNameFromHashkey('unknown'))

    def testAll(self):
        """
        The landmarkNameFromHashkey function should correctly identify the
        symbol from all landmark classes.
        """
        for cls in ALL_LANDMARK_CLASSES:
            self.assertEqual(cls.NAME, landmarkNameFromHashkey(cls.SYMBOL))

    def testProsite(self):
        """
        The landmarkNameFromHashkey function should correctly identify hashkeys
        created by the Prosite landmark class. The Prosite hashkeys are longer
        than just the class symbol.
        """
        self.assertEqual(Prosite.NAME,
                         landmarkNameFromHashkey(Prosite.SYMBOL + '00342'))


class TestAllLandmarkClasses(TestCase):
    """
    Trivial tests for the ALL_LANDMARK_CLASSES list.
    """

    def testIsAList(self):
        """
        ALL_LANDMARK_CLASSES must be a list (not a set).
        """
        self.assertTrue(isinstance(ALL_LANDMARK_CLASSES, list))

    def testAllClasses(self):
        """
        The ALL_LANDMARK_CLASSES list must be as expected.
        """
        self.assertEqual(
            [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, AminoAcids,
             BetaStrand, BetaTurn, GOR4AlphaHelix, GOR4BetaStrand, GOR4Coil,
             Prosite],
            ALL_LANDMARK_CLASSES)


class TestDefaultLandmarkClasses(TestCase):
    """
    Trivial tests for the DEFAULT_LANDMARK_CLASSES set.
    """

    def testIsAList(self):
        """
        DEFAULT_LANDMARK_CLASSES must be a list (not a set).
        """
        self.assertTrue(isinstance(DEFAULT_LANDMARK_CLASSES, list))

    def testDefaultClasses(self):
        """
        The DEFAULT_LANDMARK_CLASSES must be as expected.
        """
        self.assertEqual(
            [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, AminoAcids,
             BetaStrand, BetaTurn, Prosite],
            DEFAULT_LANDMARK_CLASSES)

    def testDefaultClassesAreInAllClasses(self):
        """
        The DEFAULT_LANDMARK_CLASSES must all appear in
        ALL_LANDMARK_CLASSES.
        """
        for klass in DEFAULT_LANDMARK_CLASSES:
            self.assertIn(klass, ALL_LANDMARK_CLASSES)
