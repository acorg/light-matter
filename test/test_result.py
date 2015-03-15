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
        self.assertIs(read, result.scannedQuery)

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
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [1],
                    'subjectLandmarkOffsets': [2],
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
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
                    'landmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'queryLandmarkOffsets': [1],
                    'queryTrigPointOffsets': [1],
                    'subjectLandmarkOffsets': [2, 10],
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
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
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [1],
                    'subjectLandmarkOffsets': [0],
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                },
                {
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [2],
                    'subjectLandmarkOffsets': [0],
                    'trigPoint': TrigPoint('Peaks', 'P', 2),
                },
                {
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [3],
                    'subjectLandmarkOffsets': [10],
                    'trigPoint': TrigPoint('Peaks', 'P', 3),
                },
            ],
        }

        result = Result(read, matches, hashCount, significanceFraction=0.25,
                        database=database)
        self.assertEqual([0], list(result.significant()))
        self.assertEqual(0.5, result.analysis[0]['bestScore'])

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
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [1],
                    'subjectLandmarkOffsets': [0],
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                },
                {
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [2],
                    'subjectLandmarkOffsets': [0],
                    'trigPoint': TrigPoint('Peaks', 'P', 2),
                },
                {
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [3],
                    'subjectLandmarkOffsets': [10],
                    'trigPoint': TrigPoint('Peaks', 'P', 3),
                },
            ],

            1: [
                {
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [1],
                    'subjectLandmarkOffsets': [0],
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                },
                {
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [2],
                    'subjectLandmarkOffsets': [0],
                    'trigPoint': TrigPoint('Peaks', 'P', 2),
                },
            ],
        }

        result = Result(read, matches, hashCount, significanceFraction=0.3,
                        database=database)
        self.assertEqual([0, 1], sorted(list(result.significant())))
        self.assertEqual(0.4, result.analysis[0]['bestScore'])
        self.assertEqual(0.4, result.analysis[1]['bestScore'])

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
        Save must write results out in the expected JSON format.
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
                    'landmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryLandmarkOffsets': [0],
                    'queryTrigPointOffsets': [1],
                    'subjectLandmarkOffsets': [0],
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                },
            ],

            1: [
                {
                    'landmark': Landmark('AlphaHelix', 'A', 27, 13),
                    'queryLandmarkOffsets': [27],
                    'queryTrigPointOffsets': [1],
                    'subjectLandmarkOffsets': [27],
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
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
                        'subjectIndex': 0,
                        'matchScore': 1.0,
                        'hsps': [
                            {
                                'hspInfo': [
                                    {
                                        'landmark': 'AlphaHelix',
                                        'landmarkLength': 9,
                                        'queryLandmarkOffset': 0,
                                        'queryTrigPointOffset': 1,
                                        'subjectLandmarkOffset': 0,
                                        'trigPoint': 'Peaks',
                                    },
                                ],
                                'score': 1.0,
                            },
                        ],
                    },
                    {
                        'subjectIndex': 1,
                        'matchScore': 1.0,
                        'hsps': [
                            {
                                'hspInfo': [
                                    {
                                        'landmark': 'AlphaHelix',
                                        'landmarkLength': 13,
                                        'queryLandmarkOffset': 27,
                                        'queryTrigPointOffset': 1,
                                        'subjectLandmarkOffset': 27,
                                        'trigPoint': 'Peaks',
                                    },
                                ],
                                'score': 1.0,
                            },
                        ],
                    },
                ],
                'queryId': 'id',
                'querySequence': 'FRRRFRRRFRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF',
            },
            result)

    def testRightNumberOfBucketsDefault(self):
        """
        If no distanceBase is specified for a database, the number of bins must
        be 20 if the length of the subject is 20 and the length of the read
        is less. This is because int(log base 1.1 of 20) = 31, which is greater
        than 20, so the value of 20 should be used. Note that 1.1 is the
        default distanceBase.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        database.addSubject(AARead('subject', 'AAAAAAAAAAAAAAAAAAAA'))
        hashCount = 1
        matches = {
            0: [
                {
                    'landmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'queryLandmarkOffsets': [1],
                    'queryTrigPointOffsets': [1],
                    'subjectLandmarkOffsets': [2],
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.0,
                        database=database, storeFullAnalysis=True)
        self.assertEqual(20, len(result.analysis[0]['histogram'].bins))

    def testRightNumberOfBucketsWithNonDefaultDistanceBase(self):
        """
        If a distanceBase of 1.3 is given and the length of the longer
        sequence (out of subject and query) is 20, there should be 11 buckets
        (because int(log base 1.3 of 20) = 11).
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [], distanceBase=1.3)
        database.addSubject(AARead('subject', 'AAAAAAAAAAAAAAAAAAAA'))
        hashCount = 1
        matches = {
            0: [
                {
                    'landmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'queryLandmarkOffsets': [1],
                    'queryTrigPointOffsets': [1],
                    'subjectLandmarkOffsets': [2],
                    'trigPoint': TrigPoint('Peaks', 'P', 1),
                },
            ],
        }
        result = Result(read, matches, hashCount, significanceFraction=0.0,
                        database=database, storeFullAnalysis=True)
        self.assertEqual(11, len(result.analysis[0]['histogram'].bins))

    def testPrintWithQueryWithNoMatchesDueToNoFinders(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to print the query and when there are no matches (in this
        case due to the database having no finders).
        """
        fp = StringIO()
        read = AARead('read', 'AGTARFSDDD')
        database = Database([], [])
        database.addSubject(read)
        result = database.find(read, significanceFraction=0.0,
                               storeFullAnalysis=True)

        result.print_(database, fp)
        expected = ("Query title: read\n"
                    "  Length: 10\n"
                    "  Covered indices: 0 (0.00%)\n"
                    "  Landmark count 0, trig point count 0\n"
                    "Overall matches: 0\n"
                    "Significant matches: 0\n"
                    "Query hash count: 0\n"
                    "Significance fraction: 0.000000\n"
                    "Significance cutoff: 0.000000\n")

        self.assertEqual(expected, fp.getvalue())

    def testPrintWithoutQueryWithNoMatchesDueToNoFinders(self):
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

        result.print_(database, fp, printQuery=False)
        expected = ("Overall matches: 0\n"
                    "Significant matches: 0\n"
                    "Query hash count: 0\n"
                    "Significance fraction: 0.000000\n"
                    "Significance cutoff: 0.000000\n")
        self.assertEqual(expected, fp.getvalue())

    def testPrintNoMatchingSubjects(self):
        """
        Check that the print_ method of a result produces the expected result
        when there are no matches.
        """
        fp = StringIO()
        read = AARead('read', 'FRRRFRRRFRFRFRFRFRFRFFRRRFRRRFRRRF')
        database = Database([AlphaHelix, BetaStrand], [])
        subject = AARead('subject', 'VICVICV')
        database.addSubject(subject)
        result = database.find(read, storeFullAnalysis=True)

        result.print_(database, fp, printQuery=False)
        expected = ("Overall matches: 0\n"
                    "Significant matches: 0\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.250000\n"
                    "Significance cutoff: 0.250000\n")
        self.assertEqual(expected, fp.getvalue())

    def testPrintOneMatchStoredAnalysis(self):
        """
        Check that the print_ method of a result produces the expected result
        when there is a match.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        read = AARead('read', sequence)
        result = database.find(read, significanceFraction=0.0,
                               storeFullAnalysis=True)

        result.print_(database, fp, printQuery=False)

        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.000000\n"
                    "Significance cutoff: 0.000000\n"
                    "Matched subjects:\n"
                    "  Subject 1 title: subject\n"
                    "    Index in database: 0\n"
                    "    Best HSP score: 1.0\n"
                    "    HSP count: 1\n"
                    "      HSP 1 has 1 matching (landmark, trigpoint) pair "
                    "and score: 1.000000\n"
                    "        Landmark AlphaHelix symbol='A' offset=0 len=9 "
                    "detail=2 subjectOffset=0\n"
                    "        Trig point AlphaHelix symbol='A' offset=25 "
                    "len=13 detail=3\n")

        self.assertEqual(expected, fp.getvalue())

    def testOneMatchNoFullAnalysis(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the read, when there are matches and the
        full analysis is not stored, and when sequences are not printed.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        read = AARead('read', sequence)
        result = database.find(read, significanceFraction=0.0)

        result.print_(database, fp, printQuery=False)
        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.000000\n"
                    "Significance cutoff: 0.000000\n"
                    "Matched subjects:\n"
                    "  Subject 1 title: subject\n"
                    "    Index in database: 0\n"
                    "    Best HSP score: 1.0\n"
                    "    HSP count: 1\n"
                    "      HSP 1 has 1 matching (landmark, trigpoint) pair "
                    "and score: 1.000000\n"
                    "        Landmark AlphaHelix symbol='A' offset=0 len=9 "
                    "detail=2 subjectOffset=0\n"
                    "        Trig point AlphaHelix symbol='A' offset=25 "
                    "len=13 detail=3\n")

        self.assertEqual(expected, fp.getvalue())

    def testOneMatchNoFullAnalysisNoFeatures(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the read, when there are matches and the
        full analysis is not stored, and when features are not printed.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        read = AARead('read', sequence)
        result = database.find(read, significanceFraction=0.0)

        result.print_(database, fp, printQuery=False, printFeatures=False)
        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.000000\n"
                    "Significance cutoff: 0.000000\n"
                    "Matched subjects:\n"
                    "  Subject 1 title: subject\n"
                    "    Index in database: 0\n"
                    "    Best HSP score: 1.0\n"
                    "    HSP count: 1\n"
                    "      HSP 1 has 1 matching (landmark, trigpoint) pair "
                    "and score: 1.000000\n")

        self.assertEqual(expected, fp.getvalue())

    def testPrintWithoutQueryWithMatchesWithoutFullAnalysisWithSequences(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the read, when there are matches and the
        full analysis is not stored, and when sequences are printed.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        read = AARead('read', sequence)
        result = database.find(read, significanceFraction=0.0)

        result.print_(database, fp, printQuery=False, printSequences=True,
                      printFeatures=False)
        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.000000\n"
                    "Significance cutoff: 0.000000\n"
                    "Matched subjects:\n"
                    "  Subject 1 title: subject\n"
                    "    Index in database: 0\n"
                    "    Best HSP score: 1.0\n"
                    "    Sequence: FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF\n"
                    "    HSP count: 1\n"
                    "      HSP 1 has 1 matching (landmark, trigpoint) pair "
                    "and score: 1.000000\n")

        self.assertEqual(expected, fp.getvalue())
