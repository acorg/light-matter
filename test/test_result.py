from unittest import TestCase

from light.result import ScannedReadDatabaseResult
from dark.reads import AARead


class TestResult(TestCase):
    """
    Tests for the light.result.ScannedReadDatabaseResult class.
    """
    def testEvaluateNotSignificantIdenticalReads(self):
        """
        A not significant result must not be returned if the matches are from
        the same reads.
        """
        read = AARead('read', 'AGTARFSDDD')
        result = ScannedReadDatabaseResult(read)
        result.addMatch({'subjectOffset': 3, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.finalize()
        self.assertEqual([], result.significant)

    def testEvaluateOneSignificant(self):
        """
        One significant result must be returned.
        """
        read = AARead('read', 'AGTARFSDDD')
        result = ScannedReadDatabaseResult(read)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 8, 'readOffset': 1}, 0, 300)
        result.finalize()
        self.assertEqual([(0, 'read', 5)], result.significant)

    def testEvaluateTwoSignificantDifferentSubjects(self):
        """
        Two significant results must be returned, when they are from different
        subjects.
        """
        read = AARead('read', 'AGTARFSDDD')
        result = ScannedReadDatabaseResult(read)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 3, 'readOffset': 1}, 0, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 1, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 1, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 1, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 1, 300)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 1, 300)
        result.addMatch({'subjectOffset': 3, 'readOffset': 1}, 1, 300)
        result.finalize()
        self.assertEqual([(0, 'read', 5), (1, 'read', 5)],
                         result.significant)
