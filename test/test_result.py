from unittest import TestCase
from cStringIO import StringIO
from json import loads

from dark.reads import AARead

from light.result import Result


class TestResult(TestCase):
    """
    Tests for the light.result.Result class.
    """
    def testNoMatches(self):
        """
        A result with no matches added to it must have the expected attributes.
        """
        read = AARead('read', 'AGTARFSDDD')
        hashCount = 0
        result = Result(read, {}, hashCount, significanceFraction=0,
                        bucketFactor=1)
        self.assertEqual({}, result.matches)
        self.assertEqual([], list(result.significant()))
        self.assertIs(read, result.read)

    def testAddOneMatch(self):
        """
        Adding information about one match should result in that information
        being stored in the result instance.
        """
        read = AARead('read', 'AGTARFSDDD')
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [2],
                    'readOffset': 1,
                    'subjectLength': 100,
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0,
                        bucketFactor=1)
        self.assertEqual(matches, result.matches)

    def testNoSignificantMatches(self):
        """
        No matches are significant if there are not enough distance deltas
        above the mean number of deltas.
        """
        read = AARead('read', 'AGTARFSDDD')
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [2, 10],
                    'readOffset': 1,
                    'subjectLength': 100,
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=5,
                        bucketFactor=1)
        self.assertEqual([], list(result.significant()))

    def testOneSignificantMatch(self):
        """
        The index of a significant result must be set correctly, and its score
        (the maximum number of identical distances) must be too.
        """
        read = AARead('read', 'AGTARFSDDD')
        hashCount = 3
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [0],
                    'readOffset': 0,
                    'subjectLength': 100,
                },
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [0],
                    'readOffset': 0,
                    'subjectLength': 100,
                },
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [10],
                    'readOffset': 0,
                    'subjectLength': 100,
                },
            ],
        }

        result = Result(read, matches, hashCount, significanceFraction=0.1,
                        bucketFactor=1)
        self.assertEqual([0], list(result.significant()))
        self.assertEqual(2, result.analysis[0]['score'])

    def testTwoSignificantMatches(self):
        """
        Two significant results must be returned, when a read matches two
        different subjects, and their scores must be correct.
        """
        read = AARead('read', 'AGTARFSDDD')
        hashCount = 5
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [0],
                    'readOffset': 0,
                    'subjectLength': 100,
                },
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [0],
                    'readOffset': 0,
                    'subjectLength': 100,
                },
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [10],
                    'readOffset': 0,
                    'subjectLength': 100,
                },
            ],

            1: [
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [0],
                    'readOffset': 0,
                    'subjectLength': 1000,
                },
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [0],
                    'readOffset': 0,
                    'subjectLength': 1000,
                },
            ],
        }

        result = Result(read, matches, hashCount, significanceFraction=0.1,
                        bucketFactor=1)
        self.assertEqual([0, 1], sorted(list(result.significant())))
        self.assertEqual(2, result.analysis[0]['score'])
        self.assertEqual(2, result.analysis[1]['score'])

    def testSaveEmpty(self):
        """
        If self.matches is empty, return an empty output.
        """
        read = AARead('read', 'AGTARFSDDD')
        result = Result(read, [], 0, 0, bucketFactor=1)
        fp = StringIO()
        result.save(fp=fp)
        result = loads(fp.getvalue())
        self.assertEqual(
            {
                'alignments': [],
                'queryId': 'read',
                'querySequence': 'AGTARFSDDD',
            },
            result)

    def testSaveReturnsItsArgument(self):
        """
        The save function must return its (fp) argument.
        """
        result = Result(AARead('id', 'A'), {}, 0, significanceFraction=0,
                        bucketFactor=1)
        fp = StringIO()
        self.assertIs(fp, result.save(fp))

    def testSave(self):
        """
        Save must produce the right JSON format.
        """
        read = AARead('id', 'FRRRFRRRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF')
        hashCount = 2
        matches = {
            0: [
                {
                    'trigPointName': 'AlphaHelix',
                    'landmarkLength': 9,
                    'readOffset': 0,
                    'subjectOffsets': [0],
                    'landmarkName': 'AlphaHelix',
                    'subjectLength': 1000,
                },
            ],

            27: [
                {
                    'trigPointName': 'AlphaHelix',
                    'landmarkLength': 13,
                    'readOffset': 27,
                    'subjectOffsets': [27],
                    'landmarkName': 'AlphaHelix',
                    'subjectLength': 100,
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.1,
                        bucketFactor=1)
        fp = StringIO()
        result.save(fp=fp)
        result = loads(fp.getvalue())
        self.assertEqual(
            {
                'alignments': [
                    {
                        'matchInfo': [
                            {
                                'landmarkLength': 9,
                                'landmarkName': 'AlphaHelix',
                                'readOffset': 0,
                                'subjectOffsets': [0],
                                'trigPointName': 'AlphaHelix',
                                'subjectLength': 1000,
                            },
                        ],
                        'matchScore': 1,
                        'subjectIndex': 0
                    },
                    {
                        'matchInfo': [
                            {
                                'landmarkLength': 13,
                                'landmarkName': 'AlphaHelix',
                                'readOffset': 27,
                                'subjectOffsets': [27],
                                'trigPointName': 'AlphaHelix',
                                'subjectLength': 100,
                            },
                        ],
                        'matchScore': 1,
                        'subjectIndex': 27,
                    },
                ],
                'queryId': 'id',
                'querySequence': 'FRRRFRRRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF',
            },
            result)

    def testRightNumberOfBucketsDefault(self):
        """
        If no bucket factor is given, the number of bins must be 11.
        """
        read = AARead('read', 'AGTARFSDDD')
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [2],
                    'readOffset': 1,
                    'subjectLength': 20,
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0,
                        bucketFactor=1,
                        storeAnalysis=True)
        self.assertEqual(20, len(result.analysis[0]['histogram']))

    def testRightNumberOfBuckets(self):
        """
        If a bucketFactor of 5 is given and the length of the longer sequence
        (out of subject and query) is 20, there should be 4 buckets.
        """
        read = AARead('read', 'AGTARFSDDD')
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffsets': [2],
                    'readOffset': 1,
                    'subjectLength': 20,
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0,
                        bucketFactor=5,
                        storeAnalysis=True)
        print result.analysis[0]['histogramBuckets']
        self.assertEqual(4, len(result.analysis[0]['histogram']))
