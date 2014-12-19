from unittest import TestCase
from cStringIO import StringIO
from json import loads

from light.result import Result
from light.database import Database
from dark.reads import AARead


class TestResult(TestCase):
    """
    Tests for the light.result.Result class.
    """
    def testEvaluateNotSignificantIdenticalReads(self):
        """
        A not significant result must not be returned if the matches are from
        the same reads.
        """
        database = Database([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addSubject(read)
        result = Result(read, database)
        result.addMatch({'subjectOffset': 3, 'readOffset': 1}, 0)
        result.addMatch({'subjectOffset': 2, 'readOffset': 1}, 0)
        result.finalize(15)
        self.assertEqual({}, result.significant)

    def testEvaluateOneNotSignificant(self):
        """
        No significant results must be returned unless the number of deltas in
        a bucket exceeds the mean by more than the passed aboveMeanThreshold
        value.
        """
        database = Database([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addSubject(read)
        result = Result(read, database)
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
        result.finalize(100)
        self.assertEqual({}, result.significant)

    def testEvaluateOneSignificant(self):
        """
        One significant result must be returned.
        """
        database = Database([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addSubject(read)
        result = Result(read, database)
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
        result.finalize(15)
        self.assertEqual(19, len(result.significant[0]['offsets']))

    def testEvaluateTwoSignificantDifferentSubjects(self):
        """
        Two significant results must be returned, when they are from different
        subjects.
        """
        database = Database([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addSubject(read)
        result = Result(read, database)
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
        result.finalize(15)
        self.assertEqual(19, len(result.significant[0]['offsets']))
        self.assertEqual(19, len(result.significant[1]['offsets']))

    def testSaveEmpty(self):
        """
        If self.matches is empty, return an empty output.
        """
        database = Database([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addSubject(read)
        result = Result(read, database)
        fp = StringIO()
        result.save(fp=fp)
        result = loads(fp.getvalue())
        self.assertEqual('read', result['query'])
        self.assertEqual([], result['alignments'])

    def testSaveReturnsItsArgument(self):
        """
        The save function must return its (fp) argument.
        """
        database = Database([], [])
        result = Result(AARead('id', 'A'), database)
        fp = StringIO()
        self.assertIs(fp, result.save(fp))

    def testSave(self):
        """
        Save must produce the right JSON format.
        """
        database = Database([], [])
        read = AARead('id', 'AGTARFSDDD')
        database.addSubject(read)
        result = Result(read, database)
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
        self.assertEqual('AGTARFSDDD', result['querySequence'])
