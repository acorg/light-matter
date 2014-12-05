from unittest import TestCase

from light.result import ScannedReadDatabaseResult


class TestResult(TestCase):
    """
    Tests for the light.result.ScannedReadDatabaseResult class.
    """
    def testEvaluateNotSignificantDifferentReads(self):
        """
        A not significant result must not be returned if the matches are from
        two different reads.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query2', 2)
        result.finalize()
        self.assertEqual([], result.significant)

    def testEvaluateNotSignificantIdenticalReads(self):
        """
        A not significant result must not be returned if the matches are from
        the same reads.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch(0, 'query', 1)
        result.addMatch(0, 'query', 2)
        result.finalize()
        self.assertEqual([], result.significant)

    def testEvaluateOneSignificant(self):
        """
        One significant result must be returned.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch(0, 'query', 1)
        result.addMatch(0, 'query', 1)
        result.addMatch(0, 'query', 1)
        result.addMatch(0, 'query', 1)
        result.addMatch(0, 'query', 1)
        result.addMatch(0, 'query', 7)
        result.finalize()
        self.assertEqual([(0, 'query', 5)], result.significant)

    def testEvaluateTwoSignificant(self):
        """
        Two significant result must be returned.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 2)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 2)
        result.finalize()
        self.assertEqual([(0, 'query2', 5), (0, 'query1', 5)],
                         result.significant)

    def testEvaluateTwoSignificantOneNotSignificant(self):
        """
        Two significant results must be returned, when one non significant
        result is present.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 2)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 1)
        result.addMatch(0, 'query2', 2)
        result.addMatch(0, 'query3', 1)
        result.addMatch(0, 'query3', 2)
        result.finalize()
        self.assertEqual([(0, 'query2', 5), (0, 'query1', 5)],
                         result.significant)

    def testEvaluateTwoSignificantDifferentSubjects(self):
        """
        Two significant results must be returned, when they are from different
        subjects.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 1)
        result.addMatch(0, 'query1', 2)
        result.addMatch(1, 'query2', 1)
        result.addMatch(1, 'query2', 1)
        result.addMatch(1, 'query2', 1)
        result.addMatch(1, 'query2', 1)
        result.addMatch(1, 'query2', 1)
        result.addMatch(1, 'query2', 2)
        result.finalize()
        self.assertEqual([(0, 'query1', 5), (1, 'query2', 5)],
                         result.significant)
