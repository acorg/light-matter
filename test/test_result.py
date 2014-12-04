from unittest import TestCase

from light.result import ScannedReadDatabaseResult


class TestResult(TestCase):
    """
    Tests for the light.result.ScannedReadDatabaseResult class.
    """
    def testResultEqual(self):
        """
        Identical results must compare equal.
        """
        result1 = ScannedReadDatabaseResult('subject', 'query', 0, 'A3:A2:23')
        result2 = ScannedReadDatabaseResult('subject', 'query', 0, 'A3:A2:23')
        self.assertEqual(result1, result2)

    def testResultsDifferingSubjectsNonEqual(self):
        """
        Results with different subjects must not compare equal.
        """
        result1 = ScannedReadDatabaseResult('subject2', 'query', 0, 'A3:A2:23')
        result2 = ScannedReadDatabaseResult('subject2', 'query', 0, 'A3:A2:23')
        self.assertNotEqual(result1, result2)

    def testResultsDifferingQueriesNonEqual(self):
        """
        Results with different queries must not compare equal.
        """
        result1 = ScannedReadDatabaseResult('subject', 'query1', 0, 'A3:A2:23')
        result2 = ScannedReadDatabaseResult('subject', 'query2', 0, 'A3:A2:23')
        self.assertNotEqual(result1, result2)

    def testResultsDifferingOffsetsNonEqual(self):
        """
        Results with different offsets must not compare equal.
        """
        result1 = ScannedReadDatabaseResult('subject', 'query', 0, 'A3:A2:23')
        result2 = ScannedReadDatabaseResult('subject', 'query', 1, 'A3:A2:23')
        self.assertNotEqual(result1, result2)

    def testResultsDifferingKeysNonEqual(self):
        """
        Results with different keys must not compare equal.
        """
        result1 = ScannedReadDatabaseResult('subject', 'query', 0, 'A3:A2:23')
        result2 = ScannedReadDatabaseResult('subject', 'query', 0, 'A3:A2:2')
        self.assertNotEqual(result1, result2)
