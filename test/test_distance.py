from unittest import TestCase
from functools import wraps

from light.distance import scale as defaultScale, _pp_scale


class TestDistanceMixin(object):
    """
    Tests for the light.distance.scale function.
    """

    def testDistanceZero(self):
        """
        If the distance is 0, the scale function must return 0.
        """
        self.assertEqual(0, self.scale(0, 2.0))

    def testBaseOne(self):
        """
        If the base is 1.0, the scale function must return the value is it
        passed.
        """
        for i in range(0, 100):
            self.assertEqual(i, self.scale(i, 1.0))

    def testBaseTwo(self):
        """
        If the base is 2.0, the scale function must return the expected result.
        """
        self.assertEqual(5, self.scale(32, 2.0))

    def testNegative(self):
        """
        If the distance is negative, the scale function must return the
        expected result.
        """
        self.assertEqual(-46, self.scale(-81, 1.1))

    def testNegativeAndPositiveOpposite(self):
        """
        If the distance is negative, the scale function must return the
        negative of what it returns for the equivalent positive value. This
        must work for small and large values (small means less than 38 with a
        base of 1.1 because int(log base 1.1 38) = 38).
        """
        self.assertEqual(self.scale(32, 1.1), -1 * self.scale(-32, 1.1))
        self.assertEqual(self.scale(320, 1.1), -1 * self.scale(-320, 1.1))

    def testResultCannotBeLarger(self):
        """
        When a distance is small, the scaled distance would be larger than the
        original if we simply took the logarithm as the scaled
        distance. Think of how the graph of y = log x versus y=x looks: for
        small x the log value is higher than x. Test that the scale
        function returns the original value in this case.
        """
        # int(Log base 1.1 of 10) = 24 so we should just have 10 returned.
        self.assertEqual(10, self.scale(10, 1.1))

    def testScaleOnePointOneSmallValues(self):
        """
        If the distanceBase is 1.1, the scale function must return the expected
        result on small values (for which int(log base 1.1 value) is greater
        than value).
        """
        self.assertEqual(0, self.scale(0, 1.1))
        self.assertEqual(0, self.scale(1, 1.1))
        for i in range(2, 39):
            self.assertEqual(i, self.scale(i, 1.1))

    def testOnePointOneRange81to88(self):
        """
        With base 1.1, distances 81-88 are all mapped to 46.
        """
        for distance in range(81, 89):
            self.assertEqual(46, self.scale(distance, 1.1))

    def testOnePointOneRange200to207(self):
        """
        With base 1.1, distances 200-207 are all mapped to 55.
        """
        for distance in range(200, 208):
            self.assertEqual(55, self.scale(distance, 1.1))

    def testOnePointOneRange208to228(self):
        """
        With base 1.1, distances 208-228 are all mapped to 56.
        """
        for distance in range(208, 229):
            self.assertEqual(56, self.scale(distance, 1.1))

    def testOnePointOneRange229to250(self):
        """
        With base 1.1, distances 229-250 are all mapped to 57.
        """
        for distance in range(229, 251):
            self.assertEqual(57, self.scale(distance, 1.1))

    def testFeatureLengthScaling(self):
        """
        With base 1.35, lengths must be scaled as expected.
        """
        # This test is present just so we have an idea of how small values
        # are scaled with a base that's (currently) the default feature
        # length base.
        for (length, scaled) in ((1, 0),
                                 (2, 2),
                                 (3, 3),
                                 (4, 4),
                                 (5, 5),
                                 (6, 5),
                                 (7, 6),
                                 (8, 6),
                                 (9, 7),
                                 (10, 7),
                                 (11, 7),
                                 (12, 8),
                                 (13, 8),
                                 (14, 8),
                                 (15, 9),
                                 (16, 9),
                                 (17, 9),
                                 (18, 9),
                                 (19, 9),
                                 (20, 9),
                                 (21, 10),
                                 (22, 10),
                                 (23, 10),
                                 (24, 10),
                                 (25, 10),
                                 (26, 10),
                                 (27, 10),
                                 (28, 11),
                                 (29, 11),
                                 (30, 11),
                                 (31, 11),
                                 (32, 11),
                                 (33, 11),
                                 (34, 11),
                                 (35, 11),
                                 (36, 11),
                                 (37, 12),
                                 (38, 12)):
            self.assertEqual(scaled, self.scale(length, 1.35))


class TestDistanceDefault(TestDistanceMixin, TestCase):
    """
    Tests for the light.distance.scale function, which is the
    C extension version (as verified below).
    """

    @wraps(defaultScale)
    def scale(self, dist, base):
        return defaultScale(dist, base)

    def testIsCExtension(self):
        """
        The scale function imported by default must be the C extension.
        """
        self.assertEqual(
            'direct call to the C function of the same name',
            self.scale.__doc__)

    def testBaseZero(self):
        """
        If the base is 0.0, the C scale function returns 0.0.

        We shouldn't ever be using a base of zero, and we'll know if we do
        (when using the C scale function) because all distances will be scaled
        to zero. Also, we check for distance base <= 0 in parameters.py and
        raise a C{ValueError} there.
        """
        self.assertEqual(0.0, self.scale(27, 0.0))


class TestDistancePurePython(TestDistanceMixin, TestCase):
    """
    Tests for the light.distance._pp_scale function, the pure Python
    version of the scale function.
    """

    @wraps(_pp_scale)
    def scale(self, dist, base):
        return _pp_scale(dist, base)

    def testNotCExtension(self):
        """
        The pure Python scale function must not be the C extension.
        """
        self.assertNotEqual(
            'direct call to the C function of the same name',
            self.scale.__doc__)

    def testBaseZero(self):
        """
        If the base is 0.0, the pure Python scale function must raise
        ValueError.
        """
        self.assertRaises(ValueError, self.scale, 27, 0.0)
