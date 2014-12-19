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
        result.addMatch({
                        'trigPointName': 'Peaks',
                        'distance': -13,
                        'landmarkLength': 9,
                        'landmarkName': 'AlphaHelix',
                        'offsets': {
                            'subjectOffset': 2,
                            'readOffset': 1,
                        }
                        }, 0)
        result.addMatch({
                        'trigPointName': 'Peaks',
                        'distance': -13,
                        'landmarkLength': 9,
                        'landmarkName': 'AlphaHelix',
                        'offsets': {
                            'subjectOffset': 2,
                            'readOffset': 1,
                        }
                        }, 0)
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
        info = [
            {
                'trigPointName': 'Peaks',
                'distance': -10,
                'landmarkLength': 9,
                'landmarkName': 'AlphaHelix',
                'offsets': {
                    'subjectOffset': 0,
                    'readOffset': 0
                }
            }
        ]
        result.matches = {0: {'info': info}}
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
        info = [
            {
                'trigPointName': 'Peaks',
                'distance': -10,
                'landmarkLength': 9,
                'landmarkName': 'AlphaHelix',
                'offsets': {
                    'subjectOffset': 0,
                    'readOffset': 0
                }
            },
            {
                'trigPointName': 'Peaks',
                'distance': -10,
                'landmarkLength': 9,
                'landmarkName': 'AlphaHelix',
                'offsets': {
                    'subjectOffset': 0,
                    'readOffset': 0
                }
            },
            {
                'trigPointName': 'Peaks',
                'distance': -10,
                'landmarkLength': 9,
                'landmarkName': 'AlphaHelix',
                'offsets': {
                    'subjectOffset': 10,
                    'readOffset': 0
                }
            }
        ]
        result.matches = {0: {'info': info}}
        result.finalize(0)
        self.assertEqual(1, len(result.significant))

    def testEvaluateTwoSignificantDifferentSubjects(self):
        """
        Two significant results must be returned, when they are from different
        subjects.
        """
        database = Database([], [])
        read = AARead('read', 'AGTARFSDDD')
        database.addSubject(read)
        result = Result(read, database)
        info = [
            {
                'trigPointName': 'Peaks',
                'distance': -10,
                'landmarkLength': 9,
                'landmarkName': 'AlphaHelix',
                'offsets': {
                    'subjectOffset': 0,
                    'readOffset': 0
                }
            },
            {
                'trigPointName': 'Peaks',
                'distance': -10,
                'landmarkLength': 9,
                'landmarkName': 'AlphaHelix',
                'offsets': {
                    'subjectOffset': 0,
                    'readOffset': 0
                }
            },
            {
                'trigPointName': 'Peaks',
                'distance': -13,
                'landmarkLength': 9,
                'landmarkName': 'AlphaHelix',
                'offsets': {
                    'subjectOffset': 10,
                    'readOffset': 0
                }
            }
        ]
        result.matches = {0: {'info': info},
                          1: {'info': info}}
        result.finalize(1)
        self.assertEqual(2, len(result.significant[0]))
        self.assertEqual(2, len(result.significant[1]))

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

    def testSave(self):
        """
        Save must produce the right JSON format.
        """
        database = Database([], [])
        read = AARead('id', 'AGTARFSDDD')
        database.addSubject(read)
        result = Result(read, database)
        result.significant = {0: {
            'offsets': [
                {
                    'readOffset': 0, 'subjectOffset': 0
                },
                {
                    'readOffset': 0, 'subjectOffset': 0
                }
            ],
            'info': [
                {
                    'distance': -10,
                    'landmarkLength': 9,
                    'landmarkName': 'A2',
                    'offsets': {
                        'readOffset': 0,
                        'subjectOffset': 0
                    },
                    'trigPointName': 'P',
                },
                {
                    'distance': -13,
                    'landmarkLength': 9,
                    'landmarkName': 'A2',
                    'offsets': {
                        'readOffset': 0,
                        'subjectOffset': 0
                    },
                    'trigPointName': 'P',
                },
            ],
            'matchScore': 15,
        }
        }
        fp = StringIO()
        result.save(fp=fp)
        result = loads(fp.getvalue())
        self.assertEqual('id', result['query'])
        self.assertEqual('AGTARFSDDD', result['querySequence'])
        self.assertEqual([{
                         'subjectIndex': 0,
                         'matchScore': 15,
                         'hsps': [
                             {
                                 'trigPointName': 'P',
                                 'distance': -10,
                                 'landmarkLength': 9,
                                 'landmarkName': 'A2',
                                 'offsets': {
                                     'subjectOffset': 0,
                                     'readOffset': 0
                                 }
                             },
                             {
                                 'trigPointName': 'P',
                                 'distance': -13,
                                 'landmarkLength': 9,
                                 'landmarkName': 'A2',
                                 'offsets': {
                                     'subjectOffset': 0,
                                     'readOffset': 0
                                 }
                             }
                         ]
                         }
        ],
            result['alignments'])
