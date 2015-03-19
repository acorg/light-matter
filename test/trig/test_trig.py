from unittest import TestCase

from light.trig import (
    findTrigPoint, findTrigPoints, trigNameFromHashkey,
    ALL_TRIG_CLASSES, DEFAULT_TRIG_CLASSES, AminoAcids, IndividualPeaks,
    IndividualTroughs, Peaks, Troughs)


class TestFindTrigPoint(TestCase):
    """
    Tests for the light.trig.findTrigPoint function.
    """

    def testFindTrigPointFails(self):
        """
        The find function should return C{None} if asked to find a trig
        point class that doesn't exist.
        """
        self.assertIs(None, findTrigPoint('silly'))

    def testFindAllClasses(self):
        """
        The find function should be able to find all trig point classes by
        name.
        """
        for klass in ALL_TRIG_CLASSES:
            self.assertIs(klass, findTrigPoint(klass.NAME))


class TestFindTrigPoints(TestCase):
    """
    Tests for the light.trig.findTrigPoints function.
    """

    def testFindFails(self):
        """
        The find function should raise a ValueError if asked to find unknown
        classes.
        """
        error = '^Unknown trig point finders: x, y\.$'
        self.assertRaisesRegexp(ValueError, error, findTrigPoints, ['x', 'y'])

    def testFindNone(self):
        """
        The find function should return an empty list if passed None.
        """
        result = findTrigPoints(None)
        self.assertEqual([], result)

    def testFindKnownClasses(self):
        """
        The find function should return known classes correctly.
        """
        result = findTrigPoints(['Peaks', 'Troughs'])
        self.assertEqual(2, len(result))
        self.assertIs(Peaks, result[0])
        self.assertIs(Troughs, result[1])


class TestTrigNameFromHashkey(TestCase):
    """
    Tests for the light.trigs.trigNameFromHashkey function.
    """

    def testFail(self):
        """
        The trigNameFromHashkey function should return None if asked to
        find a hashkey that no class created.
        """
        self.assertIs(None, trigNameFromHashkey('unknown'))

    def testAll(self):
        """
        The trigNameFromHashkey function should correctly identify the
        symbol from all trig classes.
        """
        for cls in ALL_TRIG_CLASSES:
            self.assertEqual(cls.NAME, trigNameFromHashkey(cls.SYMBOL))


class TestAllTrigClasses(TestCase):
    """
    Trivial tests for the ALL_TRIG_CLASSES set.
    """

    def testIsAList(self):
        """
        ALL_TRIG_CLASSES must be a list (not a set).
        """
        self.assertTrue(isinstance(ALL_TRIG_CLASSES, list))

    def testAllClasses(self):
        """
        The ALL_TRIG_CLASSES set must be as expected.
        """
        self.assertEqual(
            [Peaks, Troughs, AminoAcids, IndividualPeaks, IndividualTroughs],
            ALL_TRIG_CLASSES)


class TestDefaultTrigClasses(TestCase):
    """
    Trivial tests for the DEFAULT_TRIG_CLASSES set.
    """

    def testIsAList(self):
        """
        DEFAULT_TRIG_CLASSES must be a list (not a set).
        """
        self.assertTrue(isinstance(DEFAULT_TRIG_CLASSES, list))

    def testDefaultClasses(self):
        """
        The DEFAULT_TRIG_CLASSES must be as expected.
        """
        self.assertEqual(
            [Peaks, Troughs, AminoAcids],
            DEFAULT_TRIG_CLASSES)

    def testDefaultClassesAreInAllClasses(self):
        """
        The DEFAULT_TRIG_CLASSES must all appear in
        ALL_TRIG_CLASSES.
        """
        for klass in DEFAULT_TRIG_CLASSES:
            self.assertIn(klass, ALL_TRIG_CLASSES)
