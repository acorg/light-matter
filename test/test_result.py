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
        result = Result(read, {}, 0)
        self.assertEqual({}, result.matches)
        self.assertEqual(set(), result.significant)
        self.assertIs(read, result.read)

    def testAddOneMatch(self):
        """
        Adding information about one match should result in that information
        being stored in the result instance.
        """
        read = AARead('read', 'AGTARFSDDD')
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'distance': -13,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 2,
                    'readOffset': 1,
                },
            ],
        }
        result = Result(read, matches, 0)
        self.assertEqual(matches, result.matches)

    def testNoSignificantMatches(self):
        """
        No matches are significant if there are not enough distance deltas
        above the mean number of deltas.
        """
        read = AARead('read', 'AGTARFSDDD')
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'distance': -13,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 2,
                    'readOffset': 1,
                },
                {
                    'trigPointName': 'Peaks',
                    'distance': -13,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 2,
                    'readOffset': 1,
                },
            ],
        }
        result = Result(read, matches, 5)
        self.assertEqual(set(), result.significant)

    def testOneSignificantMatch(self):
        """
        The index of a significant result must be set correctly, and its score
        (the maximum number of identical distances) must be too.
        """
        read = AARead('read', 'AGTARFSDDD')
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'distance': 0,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 0,
                    'readOffset': 0,
                },
                {
                    'trigPointName': 'Peaks',
                    'distance': 0,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 0,
                    'readOffset': 0,
                },
                {
                    'trigPointName': 'Peaks',
                    'distance': -10,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 10,
                    'readOffset': 0,
                },
            ],
        }

        result = Result(read, matches, aboveMeanThreshold=0.1)
        self.assertEqual({0}, result.significant)
        self.assertEqual(2, result.scores[0])

    def testTwoSignificantMatches(self):
        """
        Two significant results must be returned, when a read matches two
        different subjects, and their scores must be correct.
        """
        read = AARead('read', 'AGTARFSDDD')
        matches = {
            0: [
                {
                    'trigPointName': 'Peaks',
                    'distance': 0,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 0,
                    'readOffset': 0,
                },
                {
                    'trigPointName': 'Peaks',
                    'distance': 0,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 0,
                    'readOffset': 0,
                },
                {
                    'trigPointName': 'Peaks',
                    'distance': -10,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 10,
                    'readOffset': 0,
                },
            ],

            1: [
                {
                    'trigPointName': 'Peaks',
                    'distance': 0,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 0,
                    'readOffset': 0,
                },
                {
                    'trigPointName': 'Peaks',
                    'distance': 0,
                    'landmarkLength': 9,
                    'landmarkName': 'AlphaHelix',
                    'subjectOffset': 0,
                    'readOffset': 0,
                },
            ],
        }

        result = Result(read, matches, aboveMeanThreshold=0.1)
        self.assertEqual({0, 1}, result.significant)
        self.assertEqual(2, result.scores[0])
        self.assertEqual(2, result.scores[1])

    def testSaveEmpty(self):
        """
        If self.matches is empty, return an empty output.
        """
        read = AARead('read', 'AGTARFSDDD')
        result = Result(read, [], 0)
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
        result = Result(AARead('id', 'A'), {}, aboveMeanThreshold=0)
        fp = StringIO()
        self.assertIs(fp, result.save(fp))

    def testSave(self):
        """
        Save must produce the right JSON format.
        """
        read = AARead('id', 'FRRRFRRRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF')
        matches = {
            0: [
                {
                    'trigPointName': 'AlphaHelix',
                    'distance': -27,
                    'landmarkLength': 9,
                    'readOffset': 0,
                    'subjectOffset': 0,
                    'landmarkName': 'AlphaHelix',
                },
            ],

            27: [
                {
                    'trigPointName': 'AlphaHelix',
                    'distance': 27,
                    'landmarkLength': 13,
                    'readOffset': 27,
                    'subjectOffset': 27,
                    'landmarkName': 'AlphaHelix',
                },
            ],
        }
        result = Result(read, matches, aboveMeanThreshold=0.1)
        fp = StringIO()
        result.save(fp=fp)
        result = loads(fp.getvalue())
        self.assertEqual(
            {
                'alignments': [
                    {
                        'matchInfo': [
                            {
                                'distance': -27,
                                'landmarkLength': 9,
                                'landmarkName': 'AlphaHelix',
                                'readOffset': 0,
                                'subjectOffset': 0,
                                'trigPointName': 'AlphaHelix'
                            },
                        ],
                        'matchScore': 1,
                        'subjectIndex': 0
                    },
                    {
                        'matchInfo': [
                            {
                                'distance': 27,
                                'landmarkLength': 13,
                                'landmarkName': 'AlphaHelix',
                                'readOffset': 27,
                                'subjectOffset': 27,
                                'trigPointName': 'AlphaHelix',
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
