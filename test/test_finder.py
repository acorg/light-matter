from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark
from light.finder import Finder
from light.landmarks import AlphaHelix, BetaStrand
from light.parameters import DatabaseParameters


class TestFinder(TestCase):
    """
    Tests for the light.finder.Finder class.
    """
    def testParamsAreStored(self):
        """
        A Finder instance must store the database parameters it is passed.
        """
        dbParams = DatabaseParameters()
        finder = Finder(dbParams)
        self.assertIs(dbParams, finder._dbParams)

    def testDefaultParams(self):
        """
        A Finder instance that is not passed any parameters must use the
        default parameters.
        """
        finder = Finder()
        self.assertIs(None, finder._dbParams.compare(DatabaseParameters()))

    def testEqual(self):
        """
        Two identical finders must compare equal.
        """
        self.assertEqual(AlphaHelix(), AlphaHelix())

    def testNotEqual(self):
        """
        Different finders must not compare equal.
        """
        self.assertNotEqual(AlphaHelix(), BetaStrand())

    def testFindWithMarginZeroAndNoMarginPresent(self):
        """
        A finder must find features when a zero margin is specified and there
        is no margin present.
        """
        read = AARead('id', 'FRRRFRRRF')
        landmark = AlphaHelix()
        result = list(landmark.findWithMargin(read, 0))
        self.assertEqual([Landmark('AlphaHelix', 'A', 0, 9, 2)], result)

    def testFindWithMarginZeroAndMarginPresentOnLeft(self):
        """
        A finder must find features when a zero margin is specified and there
        is a margin on the left.
        """
        read = AARead('id', 'MFRRRFRRRF')
        landmark = AlphaHelix()
        result = list(landmark.findWithMargin(read, 0))
        self.assertEqual([Landmark('AlphaHelix', 'A', 1, 9, 2)], result)

    def testFindWithMarginZeroAndMarginPresentOnRight(self):
        """
        A finder must find features when a zero margin is specified and there
        is a margin on the right.
        """
        read = AARead('id', 'FRRRFRRRFM')
        landmark = AlphaHelix()
        result = list(landmark.findWithMargin(read, 0))
        self.assertEqual([Landmark('AlphaHelix', 'A', 0, 9, 2)], result)

    def testFindWithMarginZeroAndMarginPresentBothSides(self):
        """
        A finder must find features when a zero margin is specified and there
        is a margin on both sides.
        """
        read = AARead('id', 'MFRRRFRRRFM')
        landmark = AlphaHelix()
        result = list(landmark.findWithMargin(read, 0))
        self.assertEqual([Landmark('AlphaHelix', 'A', 1, 9, 2)], result)

    def testFindWithMarginInsufficientMarginLeft(self):
        """
        A finder must not return a feature that has insufficient margin
        to the left.
        """
        read = AARead('id', 'MMFRRRFRRRF')
        landmark = AlphaHelix()
        result = list(landmark.findWithMargin(read, 2))
        self.assertEqual([], result)

    def testFindWithMarginInsufficientMarginRight(self):
        """
        A finder must not return a feature that has insufficient margin
        to the left.
        """
        read = AARead('id', 'FRRRFRRRFMM')
        landmark = AlphaHelix()
        result = list(landmark.findWithMargin(read, 2))
        self.assertEqual([], result)

    def testFindWithMarginSufficientMargin(self):
        """
        A finder must return a feature that has sufficient margin on both
        sides.
        """
        read = AARead('id', 'MMFRRRFRRRFMM')
        landmark = AlphaHelix()
        result = list(landmark.findWithMargin(read, 2))
        self.assertEqual([Landmark('AlphaHelix', 'A', 2, 9, 2)], result)
