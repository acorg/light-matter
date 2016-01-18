from unittest import TestCase

from dark.reads import AARead

from light.trig.random import RandomTrigPoint


class TestRandomTrigPoint(TestCase):
    """
    Tests for the TrigPoint.RandomTrigPoint class.
    """
    def testCorrectNumberOfTrigPointsDefaultDensity(self):
        """
        The correct number of TrigPoints must be returned, as specified by the
        default density parameter.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        trigPoint = RandomTrigPoint()
        result = list(trigPoint.find(read))
        self.assertEqual(2, len(result))

    def testCorrectNumberOfTrigPointsNonDefaultDensity(self):
        """
        The correct number of trigPoints must be returned, as specified by the
        density parameter.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        trigPoint = RandomTrigPoint()
        result = list(trigPoint.find(read, density=0.5))
        self.assertEqual(10, len(result))

    def testAttributes(self):
        """
        The returned trigPoints must have the correct attributes.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        trigPoint = RandomTrigPoint()
        result = list(trigPoint.find(read, density=0.05))[0]
        self.assertEqual('RandomTrigPoint', result.name)
        self.assertEqual('RT', result.symbol)
        self.assertEqual(1, result.length)
        self.assertTrue(0 <= result.offset < len(read))
