from unittest import TestCase
from cStringIO import StringIO

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

    def testSaveEmpty(self):
        """
        If self.matches is empty, return an empty output.
        """
        read = AARead('read', 'AGTARFSDDD')
        toJson = ScannedReadDatabaseResult(read)
        toJson.matches = {}
        fp = StringIO()
        toJson.save(fp=fp)
        result = fp.getvalue()
        self.assertEqual('{"query":"read","alignments":[]}\n', result)

    def testSave(self):
        """
        Save must produce the right JSON format.
        """
        read = AARead('read', 'AGTARFSDDD')
        toJson = ScannedReadDatabaseResult(read)
        toJson.matches = {0: {'absoluteOffsets':
                              [{'readOffset': 0, 'subjectOffset': 0},
                               {'readOffset': 0, 'subjectOffset': 0}],
                              'offsets': [0, 0], 'subjectLength': 15}}
        fp = StringIO()
        toJson.save(fp=fp)
        result = fp.getvalue()
        self.assertEqual('{"query":"read","alignments":[{"subjectIndex":0,'
                         '"subjectLength":15,"hsps":[{"subjectOffset":0,'
                         '"readOffset":0},{"subjectOffset":0,"readOffset":0'
                         '}]}]}\n', result)
