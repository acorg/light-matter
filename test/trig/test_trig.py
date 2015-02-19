from unittest import TestCase

from light.trig import (
    findTrigPoint, findTrigPoints, ALL_TRIG_FINDER_CLASSES,
    DEFAULT_TRIG_FINDER_CLASSES, AminoAcids, IndividualPeaks,
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
        for klass in ALL_TRIG_FINDER_CLASSES:
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


class TestAllTrigClasses(TestCase):
    """
    Trivial tests for the ALL_TRIG_FINDER_CLASSES set.
    """

    def testAllClasses(self):
        """
        The ALL_TRIG_FINDER_CLASSES set must be as expected.
        """
        self.assertEqual(
            {AminoAcids, Peaks, Troughs, IndividualPeaks, IndividualTroughs},
            ALL_TRIG_FINDER_CLASSES)


class TestDefaultTrigClasses(TestCase):
    """
    Trivial tests for the DEFAULT_TRIG_FINDER_CLASSES set.
    """

    def testDefaultClasses(self):
        """
        The DEFAULT_TRIG_FINDER_CLASSES must be as expected.
        """
        self.assertEqual(
            {Peaks, Troughs, AminoAcids},
            DEFAULT_TRIG_FINDER_CLASSES)

    def testDefaultClassesAreInAllClasses(self):
        """
        The DEFAULT_TRIG_FINDER_CLASSES must all appear in
        ALL_TRIG_FINDER_CLASSES.
        """
        for klass in DEFAULT_TRIG_FINDER_CLASSES:
            self.assertIn(klass, ALL_TRIG_FINDER_CLASSES)
