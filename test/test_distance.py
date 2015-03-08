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
        for i in range(0, 100):
            self.assertEqual(i, scale(i, 1.0))

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
        self.assertEqual(-46, scale(-81, 1.1))

    def testNegativeAndPositiveOpposite(self):
        """
        If the distance is negative, the scale function must return the
        negative of what it returns for the equivalent positive value. This
        must work for small and large values (small means less than 38 with a
        base of 1.1 because int(log base 1.1 38) = 38).
        """
        self.assertEqual(scale(32, 1.1), -1 * scale(-32, 1.1))
        self.assertEqual(scale(320, 1.1), -1 * scale(-320, 1.1))

    def testResultCannotBeLarger(self):
        """
        When a distance is small, the scaled distance would be larger than the
        original if we simply took the logarithm as the scaled
        distance. Think of how the graph of y = log x versus y=x looks: for
        small x the log value is higher than x. Test that the scale
        function returns the original value in this case.
        """
        # int(Log base 1.1 of 10) = 24 so we should just have 10 returned.
        self.assertEqual(10, scale(10, 1.1))

    def testScaleOnePointOneSmallValues(self):
        """
        If the distanceBase is 1.1, the scale function must return the expected
        result on small values (for which int(log base 1.1 value) is greater
        than value).
        """
        self.assertEqual(0, scale(0, 1.1))
        self.assertEqual(0, scale(1, 1.1))
        for i in range(2, 39):
            self.assertEqual(i, scale(i, 1.1))

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
