from unittest import TestCase

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
