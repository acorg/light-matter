from unittest import TestCase
from cStringIO import StringIO
from json import loads

from dark.reads import AARead

from light.result import Result
from light.features import Landmark, TrigPoint
from light.reads import ScannedRead
from light.database import Database
from light.landmarks import AlphaHelix, BetaStrand


class TestResult(TestCase):
    """
    Tests for the light.result.Result class.
    """
    def testNoMatches(self):
        """
        A result with no matches added to it must have the expected attributes.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        hashCount = 0
        result = Result(read, {}, hashCount, significanceFraction=0.0,
                        bucketFactor=1)
        self.assertEqual({}, result.matches)
        self.assertEqual([], list(result.significant()))
        self.assertIs(read, result.scannedRead)

    def testAddOneMatch(self):
        """
        Adding information about one match should result in that information
        being stored in the result instance.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [2],
                    'subjectLength': 100,
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.0,
                        bucketFactor=1)
        self.assertEqual(matches, result.matches)

    def testNoSignificantMatches(self):
        """
        No matches are significant if there are not enough distance deltas
        above the mean number of deltas.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'subjectOffsets': [2, 10],
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
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        hashCount = 4
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                    'subjectLength': 100,
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 2),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                    'subjectLength': 100,
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 3),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [10],
                    'subjectLength': 100,
                },
            ],
        }

        result = Result(read, matches, hashCount, significanceFraction=0.25,
                        bucketFactor=1)
        self.assertEqual([0], list(result.significant()))
        self.assertEqual(0.5, result.analysis[0]['score'])

    def testTwoSignificantMatches(self):
        """
        Two significant results must be returned, when a read matches two
        different subjects, and their scores must be correct.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        hashCount = 5
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                    'subjectLength': 100,
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 2),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                    'subjectLength': 100,
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 3),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [10],
                    'subjectLength': 100,
                },
            ],

            1: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                    'subjectLength': 1000,
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 2),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                    'subjectLength': 1000,
                },
            ],
        }

        result = Result(read, matches, hashCount, significanceFraction=0.3,
                        bucketFactor=1)
        self.assertEqual([0, 1], sorted(list(result.significant())))
        self.assertEqual(0.4, result.analysis[0]['score'])
        self.assertEqual(0.4, result.analysis[1]['score'])

    def testSaveEmpty(self):
        """
        If self.matches is empty, return an empty output.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
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
        read = ScannedRead(AARead('id', 'A'))
        result = Result(read, {}, 0, significanceFraction=0.0, bucketFactor=1)
        fp = StringIO()
        self.assertIs(fp, result.save(fp))

    def testSave(self):
        """
        Save must produce the right JSON format.
        """
        read = ScannedRead(
            AARead('id', 'FRRRFRRRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                    'subjectLength': 1000,
                },
            ],

            27: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 27, 13),
                    'subjectOffsets': [27],
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
                                'trigPointName': 'Peaks',
                                'subjectLength': 1000,
                            },
                        ],
                        'matchScore': 1.0,
                        'subjectIndex': 0
                    },
                    {
                        'matchInfo': [
                            {
                                'landmarkLength': 13,
                                'landmarkName': 'AlphaHelix',
                                'readOffset': 27,
                                'subjectOffsets': [27],
                                'trigPointName': 'Peaks',
                                'subjectLength': 100,
                            },
                        ],
                        'matchScore': 1.0,
                        'subjectIndex': 27,
                    },
                ],
                'queryId': 'id',
                'querySequence': 'FRRRFRRRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF',
            },
            result)

    def testRightNumberOfBucketsDefault(self):
        """
        If no bucket factor is given, the number of bins must be 20 if
        the length of the subject is 20 and the length of the read is less.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'subjectOffsets': [2],
                    'subjectLength': 20,
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.0,
                        bucketFactor=1, storeFullAnalysis=True)
        self.assertEqual(20, len(result.analysis[0]['histogram'].bins))

    def testRightNumberOfBuckets(self):
        """
        If a bucketFactor of 5 is given and the length of the longer sequence
        (out of subject and query) is 20, there should be 4 buckets.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'subjectOffsets': [2],
                    'subjectLength': 20,
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.0,
                        bucketFactor=5, storeFullAnalysis=True)
        self.assertEqual(4, len(result.analysis[0]['histogram'].bins))

    def testPrintWithReadWithNoMatchesDueToNoFinders(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to print the read and when there are no matches (in this
        case due to the database having no finders).
        """
        fp = StringIO()
        read = AARead('read', 'AGTARFSDDD')
        database = Database([], [])
        database.addSubject(read)
        result = database.find(read, significanceFraction=0.0,
                               storeFullAnalysis=True)

        result.print_(database, fp)
        expected = ("Read: read\n"
                    "  Sequence: AGTARFSDDD\n"
                    "  Length: 10\n"
                    "  Covered indices: 0 (0.00%)\n"
                    "  Landmark count 0, trig point count 0\n"
                    "Significant matches: 0\nOverall matches: 0\n")
        self.assertEqual(expected, fp.getvalue())

    def testPrintWithoutReadWithNoMatchesDueToNoFinders(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the read and when there are no matches (in this
        case due to the database having no finders).
        """
        fp = StringIO()
        read = AARead('read', 'AGTARFSDDD')
        database = Database([], [])
        database.addSubject(read)
        result = database.find(read, significanceFraction=0.0,
                               storeFullAnalysis=True)

        result.print_(database, fp, printRead=False)
        expected = ("Significant matches: 0\n"
                    "Overall matches: 0\n")
        self.assertEqual(expected, fp.getvalue())

    def testPrintWithoutReadWithNoMatchingSubjects(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the read and when there are no matches (in this
        case due to the database having no finders).
        """
        fp = StringIO()
        read = AARead('read', 'FRRRFRRRFRFRFRFRFRFRFFRRRFRRRFRRRF')
        database = Database([AlphaHelix, BetaStrand], [])
        subject = AARead('subject', 'VICVICV')
        database.addSubject(subject)
        result = database.find(read, storeFullAnalysis=True)

        result.print_(database, fp, printRead=False)
        expected = ("Significant matches: 0\n"
                    "Overall matches: 0\n")
        self.assertEqual(expected, fp.getvalue())

    def testPrintWithoutReadWithMatchesFullAnalysis(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the read and when there are matches and the
        full analysis is stored.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        read = AARead('read', sequence)
        result = database.find(read, significanceFraction=0.0,
                               storeFullAnalysis=True)

        result.print_(database, fp, printRead=False)
        expected = ("Significant matches: 1\n"
                    "Overall matches: 1\n"
                    "Subject matches:\n"
                    "  Title: subject\n"
                    "    Score: 1.0\n"
                    "    Sequence: FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF\n"
                    "    Database subject index: 0\n"
                    "    Hash count: 1\n"
                    "    Histogram\n"
                    "      Number of bins: 38\n"
                    "      Significance cutoff: 0.0\n"
                    "      Significant bin count: 1\n"
                    "      Max bin count: 1\n"
                    "      Counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "
                    "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "
                    "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]\n"
                    "      Max offset delta: 0\n"
                    "      Min offset delta: 0\n")

        self.assertEqual(expected, fp.getvalue())

    def testPrintWithoutReadWithMatchesWithoutFullAnalysis(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the read and when there are matches and the
        full analysis is not stored.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        read = AARead('read', sequence)
        result = database.find(read, significanceFraction=0.0)

        result.print_(database, fp, printRead=False)
        expected = ("Significant matches: 1\n"
                    "Overall matches: 1\n"
                    "Subject matches:\n"
                    "  Title: subject\n"
                    "    Score: 1.0\n"
                    "    Sequence: FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF\n"
                    "    Database subject index: 0\n")

        self.assertEqual(expected, fp.getvalue())
