from unittest import TestCase

import light.performance


class TestPerformance(TestCase):
    """
    Tests for the light.performance module (i.e., its __init__.py).
    """

    def testTestArgs(self):
        """
        The value of light.performance.testArgs must be None.
        """
        self.assertIs(None, light.performance.testArgs)
