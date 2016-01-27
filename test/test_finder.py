from unittest import TestCase

from light.features import _Feature


class TestFeature(TestCase):
    """
    Tests for the light.features._Feature class
    """
    def testZeroCoverage(self):
        """
        A feature of zero length must not cover any AA offsets.
        """
        self.assertEqual(set(),
                         _Feature('name', 'symbol', 10, 0).coveredOffsets())

    def testCoverage(self):
        """
        The coveredOffsets function must give the expected result for a
        non-zero length feature.
        """
        self.assertEqual(set([10, 11, 12, 13, 14]),
                         _Feature('name', 'symbol', 10, 5).coveredOffsets())
