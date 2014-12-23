from unittest import TestCase

from light.trig import (
    findTrigPoint, ALL_TRIG_FINDER_CLASSES, DEFAULT_TRIG_FINDER_CLASSES)
from light.trig.peaks import Peaks
from light.trig.troughs import Troughs
from light.trig.amino_acids import AminoAcids
from light.trig.individual_peaks import IndividualPeaks
from light.trig.individual_troughs import IndividualTroughs


class TestFindTrig(TestCase):
    """
    Tests for the light.trig.find function.
    """

    def testFindTrigFails(self):
        """
        The find function should return C{None} if asked to find a trig
        class that doesn't exist.
        """
        self.assertIs(None, findTrigPoint('silly'))

    def testFindPeaks(self):
        """
        The find function should be able to find all trig point classes by
        name.
        """
        for klass in ALL_TRIG_FINDER_CLASSES:
            self.assertIs(klass, findTrigPoint(klass.NAME))


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
