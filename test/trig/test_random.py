from unittest import TestCase

from dark.reads import AARead

from light.parameters import DatabaseParameters
from light.trig.random import RandomTrigPoint


class TestRandomTrigPoint(TestCase):
    """
    Tests for the TrigPoint.RandomTrigPoint class.
    """
    def testCorrectNumberOfTrigPointsDefaultDensity(self):
        """
        The correct number of trig points must be returned, based on the
        default density parameter value (0.1).
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        trigPoint = RandomTrigPoint()
        result = list(trigPoint.find(read))
        self.assertEqual(2, len(result))

    def testCorrectNumberOfTrigPointsNonDefaultDensity(self):
        """
        The correct number of trig points must be returned, as specified by the
        density parameter.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        dbParams = DatabaseParameters(randomTrigPointDensity=0.5)
        trigPoint = RandomTrigPoint(dbParams)
        result = list(trigPoint.find(read))
        self.assertEqual(10, len(result))

    def testAttributes(self):
        """
        The returned trig points must have the correct attributes.
        """
        read = AARead('id', 'FRFRFRFRFRFRFRFRFRFF')
        # Use an extremely high density, to make sure the finder never
        # fails to produce a random trig point.
        dbParams = DatabaseParameters(randomTrigPointDensity=0.9)
        trigPoint = RandomTrigPoint(dbParams)
        result = list(trigPoint.find(read))[0]
        self.assertEqual('RandomTrigPoint', result.name)
        self.assertEqual('RT', result.symbol)
        self.assertEqual(1, result.length)
        self.assertTrue(0 <= result.offset < len(read))
