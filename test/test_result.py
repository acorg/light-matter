from unittest import TestCase
from cStringIO import StringIO
from json import loads

from light.result import ScannedReadDatabaseResult
from light.database import ScannedReadDatabase
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
        database = ScannedReadDatabase([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addRead(read)
        result = ScannedReadDatabaseResult(read, database)
        result.addMatch({'subjectOffset': 3, 'readOffset': 1}, 0)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0)
        result.finalize()
        self.assertEqual({}, result.significant)

    def testEvaluateOneSignificant(self):
        """
        One significant result must be returned.
        """
        database = ScannedReadDatabase([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addRead(read)
        result = ScannedReadDatabaseResult(read, database)
        offsets = [{'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 3, 'readOffset': 1}]
        result.matches = {0: {'offsets': offsets}}
        result.finalize()
        self.assertEqual(19, len(result.significant[0]['offsets']))

    def testEvaluateTwoSignificantDifferentSubjects(self):
        """
        Two significant results must be returned, when they are from different
        subjects.
        """
        database = ScannedReadDatabase([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addRead(read)
        result = ScannedReadDatabaseResult(read, database)
        offsets = [{'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 2, 'readOffset': 1},
                   {'subjectOffset': 3, 'readOffset': 1}]
        result.matches = {0: {'offsets': offsets},
                          1: {'offsets': offsets}}
        result.finalize()
        self.assertEqual(19, len(result.significant[0]['offsets']))
        self.assertEqual(19, len(result.significant[1]['offsets']))

    def testSaveEmpty(self):
        """
        If self.matches is empty, return an empty output.
        """
        database = ScannedReadDatabase([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addRead(read)
        result = ScannedReadDatabaseResult(read, database)
        fp = StringIO()
        result.save(fp=fp)
        result = loads(fp.getvalue())
        self.assertEqual('read', result['query'])
        self.assertEqual([], result['alignments'])

    def testSave(self):
        """
        Save must produce the right JSON format.
        """
        database = ScannedReadDatabase([], [])
        read = AARead('id', 'AGTARFSDDD')
        database.addRead(read)
        result = ScannedReadDatabaseResult(read, database)
        result.significant = {0: {'offsets':
                                  [{'readOffset': 0, 'subjectOffset': 0},
                                   {'readOffset': 0, 'subjectOffset': 0}],
                                  'matchScore': 15}}
        fp = StringIO()
        result.save(fp=fp)
        result = loads(fp.getvalue())
        self.assertEqual('id', result['query'])
        self.assertEqual(
            [
                {
                    'subjectIndex': 0,
                    'subjectLength': 10,
                    'subjectTitle': 'id',
                    'hsps': [
                        {
                            'subjectOffset': 0,
                            'readOffset': 0
                        },
                        {
                            'subjectOffset': 0,
                            'readOffset': 0
                        }
                    ],
                    'matchScore': 15
                }
            ],
            result['alignments'])
