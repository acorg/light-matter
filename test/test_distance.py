from unittest import TestCase

from light.distance import scale


class TestDistance(TestCase):
    """
    Tests for the light.distance.scale function.
    """

    def testDistanceZero(self):
        """
        If the distance is 0, the scale function must return 0.
        """
        self.assertEqual(0, scale(0, 2.0))

    def testBaseZero(self):
        """
        If the base is 0.0, the scale function must raise ValueError.
        """
        self.assertRaises(ValueError, scale, 27, 0.0)

    def testBaseOne(self):
        """
        If the base is 1.0, the scale function must return the value is it
        passed.
        """
        self.assertEqual(88, scale(88, 1.0))

    def testBaseTwo(self):
        """
        If the base is 2.0, the scale function must return the expected result.
        """
        self.assertEqual(5, scale(32, 2.0))

    def testNegative(self):
        """
        If the distance is negative, the scale function must return the
        expected result.
        """
        self.assertEqual(-36, scale(-32, 1.1))

    def testNegativeAndPositiveOpposite(self):
        """
        If the distance is negative, the scale function must return the
        negative of what it returns for the equivalent positive value.
        """
        self.assertEqual(scale(32, 1.1), -1 * scale(-32, 1.1))

    def testResultCanBeLarger(self):
        """
        When a distance is small, the scaled distance can be larger. This may
        seem a little weird, but it's harmless. The scaled values quickly
        become less than the passed values. Think of how the graph of y = log x
        versus y=x looks: for small x the log value is higher than x.
        """
        self.assertEqual(24, scale(10, 1.1))

    def testScaleOnePointOne(self):
        """
        If the distanceBase is 1.1, the scale function must return the expected
        values.
        """
        for distance, expected in (
                (0, 0), (1, 0), (2, 7), (3, 11), (4, 14), (5, 16), (6, 18),
                (7, 20), (8, 21), (9, 23), (10, 24), (11, 25), (12, 26),
                (13, 26), (14, 27)):
            self.assertEqual(expected, scale(distance, 1.1))

    def testOnePointOneRange81to88(self):
        """
        With base is 1.1, distances 81-88 are all mapped to 46.
        """
        for distance in range(81, 89):
            self.assertEqual(46, scale(distance, 1.1))

    def testOnePointOneRange200to207(self):
        """
        With base is 1.1, distances 200-207 are all mapped to 55.
        """
        for distance in range(200, 208):
            self.assertEqual(55, scale(distance, 1.1))

    def testOnePointOneRange208to228(self):
        """
        With base is 1.1, distances 208-228 are all mapped to 56.
        """
        for distance in range(208, 229):
            self.assertEqual(56, scale(distance, 1.1))

    def testOnePointOneRange229to250(self):
        """
        With base is 1.1, distances 229-250 are all mapped to 57.
        """
        for distance in range(229, 251):
            self.assertEqual(57, scale(distance, 1.1))
