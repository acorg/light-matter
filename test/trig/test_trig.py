from unittest import TestCase

from light.trig import find, ALL_TRIG_FINDER_CLASSES
from light.trig.peaks import Peaks
from light.trig.troughs import Troughs


class TestFindTrig(TestCase):
    """
    Tests for the light.trig.find function.
    """

    def testFindTrigFails(self):
        """
        The find function should return C{None} if asked to find a trig
        class that doesn't exist.
        """
        self.assertIs(None, find('silly'))

    def testFindPeaks(self):
        """
        The find function should be able to find the Peak class by name.
        """
        self.assertIs(Peaks, find('Peaks'))


class TestAllTrigClasses(TestCase):
    """
    Trivial tests for the ALL_TRIG_FINDER_CLASSES set.
    """

    def testAllClassesContainsTroughs(self):
        """
        The ALL_TRIG_FINDER_CLASSES set must contain the Troughs class.
        """
        self.assertIn(Troughs, ALL_TRIG_FINDER_CLASSES)

    def testAllClassesContainsPeaks(self):
        """
        The ALL_LANDMARK_FINDER_CLASSES set must contain the Peaks
        class.
        """
        self.assertIn(Peaks, ALL_TRIG_FINDER_CLASSES)
