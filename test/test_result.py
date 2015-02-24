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
        database = Database([], [])
        hashCount = 0
        result = Result(read, {}, hashCount, significanceFraction=0.0,
                        database=database)
        self.assertEqual({}, result.matches)
        self.assertEqual([], list(result.significant()))
        self.assertIs(read, result.scannedRead)

    def testAddOneMatch(self):
        """
        Adding information about one match should result in that information
        being stored in the result instance.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        database.addSubject(AARead('subject', 'AAA'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [2],
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.0,
                        database=database)
        self.assertEqual(matches, result.matches)

    def testNoSignificantMatches(self):
        """
        No matches are significant if there are not enough distance deltas
        above the mean number of deltas.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        database.addSubject(AARead('subject', 'AAA'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'subjectOffsets': [2, 10],
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=5,
                        database=database)
        self.assertEqual([], list(result.significant()))

    def testOneSignificantMatch(self):
        """
        The index of a significant result must be set correctly, and its score
        (the maximum number of identical distances) must be too.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        database.addSubject(AARead('subject', 'AAA'))
        hashCount = 4
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 2),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 3),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [10],
                },
            ],
        }

        result = Result(read, matches, hashCount, significanceFraction=0.25,
                        database=database)
        self.assertEqual([0], list(result.significant()))
        self.assertEqual(0.5, result.analysis[0]['score'])

    def testTwoSignificantMatches(self):
        """
        Two significant results must be returned, when a read matches two
        different subjects, and their scores must be correct.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        database.addSubject(AARead('subject0', 'AAA'))
        database.addSubject(AARead('subject1', 'AAA'))
        hashCount = 5
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 2),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 3),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [10],
                },
            ],

            1: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                },
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 2),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                },
            ],
        }

        result = Result(read, matches, hashCount, significanceFraction=0.3,
                        database=database)
        self.assertEqual([0, 1], sorted(list(result.significant())))
        self.assertEqual(0.4, result.analysis[0]['score'])
        self.assertEqual(0.4, result.analysis[1]['score'])

    def testSaveEmpty(self):
        """
        If self.matches is empty, return an empty output.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        result = Result(read, [], 0, 0, database=database)
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
        database = Database([], [])
        result = Result(read, {}, 0, significanceFraction=0.0,
                        database=database)
        fp = StringIO()
        self.assertIs(fp, result.save(fp))

    def testSave(self):
        """
        Save must produce the right JSON format.
        """
        read = ScannedRead(
            AARead('id', 'FRRRFRRRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'))
        database = Database([], [])
        database.addSubject(AARead('subject0', 'AAA'))
        database.addSubject(AARead('subject1', 'AAA'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectOffsets': [0],
                },
            ],

            1: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 27, 13),
                    'subjectOffsets': [27],
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.1,
                        database=database)
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
                            },
                        ],
                        'matchScore': 1.0,
                        'subjectIndex': 1,
                    },
                ],
                'queryId': 'id',
                'querySequence': 'FRRRFRRRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF',
            },
            result)

    def testRightNumberOfBucketsDefault(self):
        """
        If no distanceBase is specified for a database, the number of bins must
        be 31 if the length of the subject is 20 and the length of the read
        is less. This is because log base 1.1 of 20 is 31.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        database.addSubject(AARead('subject', 'AAAAAAAAAAAAAAAAAAAA'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'subjectOffsets': [2],
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.0,
                        database=database, storeFullAnalysis=True)
        self.assertEqual(31, len(result.analysis[0]['histogram'].bins))

    def testRightNumberOfBucketsWithNonDefaultDistanceBase(self):
        """
        If a distanceBase of 1.3 is given and the length of the longer
        sequence (out of subject and query) is 20, there should be 11 buckets
        (because int(log base 1.3 of 20) = 11).
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        database.addSubject(AARead('subject', 'AAAAAAAAAAAAAAAAAAAA'))
        hashCount = 1
        matches = {
            0: [
                {
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                    'landmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'subjectOffsets': [2],
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.0,
                        database=database, storeFullAnalysis=True)
        self.assertEqual(31, len(result.analysis[0]['histogram'].bins))

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
