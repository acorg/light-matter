import six
from unittest import TestCase

from light.performance.stats import Stats


class TestStats(TestCase):
    """
    Tests for the light.performance.stats.Stats class.
    """
    def testNoValuesSummary(self):
        """
        If no values have been added to a Stats instance, its summary
        method must raise ValueError.
        """
        s = Stats()
        error = '^No values have been added to Stats instance\.$'
        six.assertRaisesRegex(self, ValueError, error, s.summary)

    def testNoValuesString(self):
        """
        If no values have been added to a Stats instance, its __str__
        method must raise ValueError.
        """
        s = Stats()
        error = '^No values have been added to Stats instance\.$'
        six.assertRaisesRegex(self, ValueError, error, str, s)

    def testEqual(self):
        """
        Stats instances that have had the same values added to them in
        the same order must be considered equal.
        """
        s1 = Stats()
        s1.extend([3, 4, 8])
        s2 = Stats()
        s2.extend([3, 4, 8])
        self.assertEqual(s1, s2)

    def testUnequal(self):
        """
        Stats instances that have had the same values added to them, but in
        a different order must be considered unequal.
        """
        s1 = Stats()
        s1.extend([3, 4, 8])
        s2 = Stats()
        s2.extend([3, 8, 4])
        self.assertNotEqual(s1, s2)

    def testOneValue(self):
        """
        If one value has been added to a Stats instance, its summary
        method must return the expected result.
        """
        s = Stats()
        s.append(3)
        self.assertEqual(
            {
                'count': 1,
                'max': 3,
                'mean': 3.0,
                'median': 3.0,
                'min': 3,
                'sd': 0.0,
                'sum': 3,
                'variance': 0.0,
            },
            s.summary())

    def testTwoValues(self):
        """
        If two values have been added to a Stats instance, its summary
        method must return the expected result.
        """
        s = Stats()
        s.append(3)
        s.append(4)
        self.assertEqual(
            {
                'count': 2,
                'max': 4,
                'mean': 3.5,
                'median': 3.5,
                'min': 3,
                'sd': 0.5,
                'sum': 7,
                'variance': 0.25,
            },
            s.summary())

    def testTwoValuesString(self):
        """
        If two values have been added to a Stats instance, its __str__
        method must return the expected result.
        """
        s = Stats()
        s.append(3)
        s.append(4)
        self.assertEqual(('Summary:\n'
                          '  Count: 2\n'
                          '  Max: 4\n'
                          '  Mean: 3.5000\n'
                          '  Median: 3.5000\n'
                          '  Min: 3\n'
                          '  SD: 0.5000\n'
                          '  Variance: 0.2500'),
                         str(s))

    def testThreeValues(self):
        """
        If three values have been added to a Stats instance, its summary
        method must return the expected result.
        """
        s = Stats()
        s.append(3)
        s.append(4)
        s.append(8)
        result = s.summary()

        self.assertAlmostEqual(2.16024689, result['sd'])
        self.assertAlmostEqual(4.66666666, result['variance'])
        del result['sd']
        del result['variance']

        self.assertEqual(
            {
                'count': 3,
                'max': 8,
                'mean': 5.0,
                'median': 4.0,
                'min': 3,
                'sum': 15,
            },
            result)

    def testThreeValuesString(self):
        """
        If three values have been added to a Stats instance, its __str__
        method must return the expected result.
        """
        s = Stats()
        s.append(3)
        s.append(4)
        s.append(8)
        self.assertEqual(('Summary:\n'
                          '  Count: 3\n'
                          '  Max: 8\n'
                          '  Mean: 5.0000\n'
                          '  Median: 4.0000\n'
                          '  Min: 3\n'
                          '  SD: 2.1602\n'
                          '  Variance: 4.6667'),
                         str(s))

    def testThreeValuesViaExtend(self):
        """
        If three values have been added to a Stats instance via extend, its
        summary method must return the expected result.
        """
        s = Stats()
        s.extend([3, 4, 8])
        result = s.summary()

        self.assertAlmostEqual(2.16024689, result['sd'])
        self.assertAlmostEqual(4.66666666, result['variance'])
        del result['sd']
        del result['variance']

        self.assertEqual(
            {
                'count': 3,
                'max': 8,
                'mean': 5.0,
                'median': 4.0,
                'min': 3,
                'sum': 15,
            },
            result)

    def testDescription(self):
        """
        If a description string is given to Stats.__init__, it must appear
        in the __str__ result.
        """
        s = Stats('Arm length')
        s.append(3)
        s.append(4)
        s.append(8)
        self.assertEqual(('Arm length:\n'
                          '  Count: 3\n'
                          '  Max: 8\n'
                          '  Mean: 5.0000\n'
                          '  Median: 4.0000\n'
                          '  Min: 3\n'
                          '  SD: 2.1602\n'
                          '  Variance: 4.6667'),
                         str(s))
