from unittest import TestCase

from light.result import ScannedReadDatabaseResult


class TestResult(TestCase):
    """
    Tests for the light.result.ScannedReadDatabaseResult class.
    """
    def testEvaluateNotSignificantDifferentReads(self):
        """
        A not significant result must not be returned.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query2', 2)
        result.finalize()
        self.assertEqual([], result.significant)

    def testEvaluateNotSignificantIdenticalReads(self):
        """
        A not significant result must not be returned.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch('subject', 'query', 1)
        result.addMatch('subject', 'query', 2)
        result.finalize()
        self.assertEqual([], result.significant)

    def testEvaluateOneSignificant(self):
        """
        One significant result must be returned.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch('subject', 'query', 1)
        result.addMatch('subject', 'query', 1)
        result.addMatch('subject', 'query', 1)
        result.addMatch('subject', 'query', 1)
        result.addMatch('subject', 'query', 1)
        result.addMatch('subject', 'query', 7)
        result.finalize()
        self.assertEqual([('subject', 'query', 5)], result.significant)

    def testEvaluateTwoSignificant(self):
        """
        Two significant result must be returned.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 2)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 2)
        result.finalize()
        self.assertEqual([('subject', 'query2', 5), ('subject', 'query1', 5)],
                         result.significant)

    def testEvaluateTwoSignificantOneNotSignificant(self):
        """
        Two significant result must be returned, when one non-significant
        result is present.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 1)
        result.addMatch('subject', 'query1', 2)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 1)
        result.addMatch('subject', 'query2', 2)
        result.addMatch('subject', 'query3', 1)
        result.addMatch('subject', 'query3', 2)
        result.finalize()
        self.assertEqual([('subject', 'query2', 5), ('subject', 'query1', 5)],
                         result.significant)

    def testEvaluateTwoSignificantDifferentSubjects(self):
        """
        Two significant result must be returned, when they are from different
        subjects.
        """
        result = ScannedReadDatabaseResult()
        result.addMatch('subject1', 'query1', 1)
        result.addMatch('subject1', 'query1', 1)
        result.addMatch('subject1', 'query1', 1)
        result.addMatch('subject1', 'query1', 1)
        result.addMatch('subject1', 'query1', 1)
        result.addMatch('subject1', 'query1', 2)
        result.addMatch('subject2', 'query2', 1)
        result.addMatch('subject2', 'query2', 1)
        result.addMatch('subject2', 'query2', 1)
        result.addMatch('subject2', 'query2', 1)
        result.addMatch('subject2', 'query2', 1)
        result.addMatch('subject2', 'query2', 2)
        result.finalize()
        self.assertEqual([('subject1', 'query1', 5),
                          ('subject2', 'query2', 5)], result.significant)
