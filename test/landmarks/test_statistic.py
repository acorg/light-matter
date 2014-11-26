from unittest import TestCase

from light.landmarks import find
from light.landmarks.alpha_helix import AlphaHelix


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
