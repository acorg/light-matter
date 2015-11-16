from unittest import TestCase
from io import StringIO
from json import loads
import warnings

from dark.reads import AARead

from light.result import Result
from light.parameters import Parameters, FindParameters
from light.features import Landmark, TrigPoint
from light.database import Database
from light.landmarks import AlphaHelix, BetaStrand, AminoAcids as AminoAcidsLm
from light.trig import Troughs, AminoAcids


class TestResult(TestCase):
    """
    Tests for the light.result.Result class.
    """
    def testNoMatches(self):
        """
        A result with no matches added to it must have the expected attributes.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [])
        database = Database(params)
        hashCount = 0
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.1)
        result = Result(read, database, {}, hashCount, findParams)
        self.assertEqual({}, result.matches)
        self.assertEqual([], list(result.significantSubjects()))
        self.assertIs(read, result.query)

    def testUnknownSignificanceMethod(self):
        """
        The Result class must raise ValueError if passed an unknown
        significance method.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [])
        database = Database(params)
        database.addSubject(AARead('subject', 'AAA'))
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 1),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 1),
                },
            ],
        }
        hashCount = 1
        error = "^Unknown significance method 'xxx'$"
        findParams = FindParameters(significanceMethod='xxx',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.1)
        self.assertRaisesRegex(
            ValueError, error, Result, read, database._connector,
            matches, hashCount, findParams)

    def testAddOneMatch(self):
        """
        Passing information in which one match is present should result in
        that information being stored in the result instance.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [])
        database = Database(params)
        database.addSubject(AARead('subject', 'AAA'))
        hashCount = 1
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 1),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 1),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.1)
        result = Result(read, database, matches, hashCount, findParams)
        self.assertEqual(matches, result.matches)

    def testNoSignificantMatches(self):
        """
        No matches are significant if there are not enough distance deltas
        above the mean number of deltas.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([AlphaHelix], [AminoAcids])
        database = Database(params)
        database.addSubject(AARead('subject', 'ADDDADDDAWW'))
        hashCount = 1
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 1),
                    'queryTrigPoint': TrigPoint('AminoAcids', 'M', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 8),
                    'subjectTrigPoint': TrigPoint('AminoAcids', 'M', 9),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 1),
                    'queryTrigPoint': TrigPoint('AminoAcids', 'M', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('AminoAcids', 'M', 10),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=5)
        result = Result(read, database, matches, hashCount, findParams)
        self.assertEqual([], list(result.significantSubjects()))

    def testOneSignificantMatch(self):
        """
        The index of a significant result must be set correctly, and its score
        (the maximum number of identical distances) must be too.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([AlphaHelix], [AminoAcids])
        database = Database(params)
        database.addSubject(AARead('subject', 'ADDDADDDAWWWW'))
        hashCount = 4
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 2),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 3),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 10, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.25)
        result = Result(read, database, matches, hashCount, findParams)
        self.assertEqual(['0'], list(result.significantSubjects()))
        self.assertEqual(0.5, result.analysis['0']['bestScore'])

    def testTwoSignificantMatches(self):
        """
        Two significant results must be returned, when a read matches two
        different subjects, and their scores must be correct.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([AlphaHelix], [AminoAcids])
        database = Database(params)
        # Note that both subject added here have a hash count of 5 (the
        # same as the passed query hash count). If we use anything less,
        # the score will change because its calculation uses the min
        # query/subject hash count in the denominator.
        database.addSubject(AARead('subject0', 'ADDDADDDAWWWWW'))
        database.addSubject(AARead('subject1', 'ADDDADDDAWWWWW'))
        hashCount = 5
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 0),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 1),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 2),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 3),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 10, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],

            '1': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 2),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.3)
        result = Result(read, database, matches, hashCount, findParams)
        self.assertEqual(['0', '1'],
                         sorted(list(result.significantSubjects())))
        self.assertEqual(0.4, result.analysis['0']['bestScore'])
        self.assertEqual(0.4, result.analysis['1']['bestScore'])

    def testSaveEmpty(self):
        """
        If self.matches is empty, return an empty output.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [])
        database = Database(params)
        result = Result(read, database, [], 0, 0, 'HashFraction',
                        'MinHashesScore')
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
        read = AARead('id', 'A')
        params = Parameters([], [])
        database = Database(params)
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.1)
        result = Result(read, database, {}, 0, findParams)
        fp = StringIO()
        self.assertIs(fp, result.save(fp))

    def testSave(self):
        """
        Save must write results out in the expected JSON format.
        """
        read = AARead('id', 'FRRRFRRRFRFRFRFRRRFRRRF')
        params = Parameters([AlphaHelix], [])
        database = Database(params)
        database.addSubject(AARead('subject0', 'FRRRFRRRFRFRFRFRRRFRRRF'))
        database.addSubject(AARead('subject1', 'FRRRFRRRFRFRFRFRRRFRRRF'))
        hashCount = 1
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],

            '1': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 27, 13),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 27, 13),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.1)
        result = Result(read, database, matches, hashCount, findParams)
        fp = StringIO()
        result.save(fp=fp)
        result = loads(fp.getvalue())
        self.assertEqual(
            {
                'alignments': [
                    {
                        'subjectIndex': '0',
                        'matchScore': 1.0,
                        'hsps': [
                            {
                                'hspInfo': [
                                    {
                                        'landmark': 'AlphaHelix',
                                        'queryLandmarkLength': 9,
                                        'queryLandmarkOffset': 0,
                                        'subjectLandmarkLength': 9,
                                        'subjectLandmarkOffset': 0,
                                        'trigPoint': 'Peaks',
                                        'queryTrigPointOffset': 1,
                                        'subjectTrigPointOffset': 0,
                                    },
                                ],
                                'score': 1.0,
                            },
                        ],
                    },
                    {
                        'subjectIndex': '1',
                        'matchScore': 1.0,
                        'hsps': [
                            {
                                'hspInfo': [
                                    {
                                        'landmark': 'AlphaHelix',
                                        'queryLandmarkLength': 13,
                                        'queryLandmarkOffset': 27,
                                        'subjectLandmarkLength': 13,
                                        'subjectLandmarkOffset': 27,
                                        'trigPoint': 'Peaks',
                                        'queryTrigPointOffset': 1,
                                        'subjectTrigPointOffset': 0,
                                    },
                                ],
                                'score': 1.0,
                            },
                        ],
                    },
                ],
                'queryId': 'id',
                'querySequence': 'FRRRFRRRFRFRFRFRRRFRRRF',
            },
            result)

    def testRightNumberOfBucketsDefault(self):
        """
        If no distanceBase is specified for a database, the number of bins must
        be 21 if the length of the subject is 21 and the length of the read
        is less. This is because int(log base 1.1 of 21) = 31, which is
        greater than 21, so the value of 21 should be used. Note that 1.1
        is the default distanceBase.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [])
        database = Database(params)
        database.addSubject(AARead('subject', 'A' * 21))
        hashCount = 1
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 2, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.1)
        result = Result(read, database, matches, hashCount, findParams,
                        storeFullAnalysis=True)
        self.assertEqual(21, result.analysis['0']['histogram'].nBins)

    def testRightNumberOfBucketsDefaultNonEven(self):
        """
        If no distanceBase is specified for a database, the number of bins must
        be 21 if the length of the subject is 20 and the length of the read
        is less. This is because int(log base 1.1 of 20) = 31, which is greater
        than 20, so the value of 20 should be used but result.py will adjust
        this to 21 to ensure that the number of bins is odd.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [])
        database = Database(params)
        database.addSubject(AARead('subject', 'A' * 20))
        hashCount = 1
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 2, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.1)
        result = Result(read, database, matches, hashCount, findParams,
                        storeFullAnalysis=True)
        self.assertEqual(21, result.analysis['0']['histogram'].nBins)

    def testRightNumberOfBucketsWithNonDefaultDistanceBase(self):
        """
        If a distanceBase of 1.3 is given and the length of the longer
        sequence (out of subject and query) is 20, there should be 11 buckets
        (because int(log base 1.3 of 20) = 11).
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [], distanceBase=1.3)
        database = Database(params)
        database.addSubject(AARead('subject', 'AAAAAAAAAAAAAAAAAAAA'))
        hashCount = 1
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 2, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.1)
        result = Result(read, database, matches, hashCount, findParams,
                        storeFullAnalysis=True)
        self.assertEqual(11, result.analysis['0']['histogram'].nBins)

    def testPrintWithQueryWithNoMatchesDueToNoFinders(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to print the query and when there are no matches (in this
        case due to the database having no finders).
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [])
        database = Database(params)
        database.addSubject(read)
        findParams = FindParameters(significanceFraction=0.1,
                                    scoreMethod='MinHashesScore')
        result = database.find(read, findParams, storeFullAnalysis=True)

        expected = ('Query title: read\n'
                    '  Length: 10\n'
                    '  Covered indices: 0 (0.00%)\n'
                    '  Landmark count 0, trig point count 0\n'
                    'Find parameters:\n'
                    '  Significance method: HashFraction\n'
                    '  Significance fraction: 0.100000\n'
                    '  Score method: MinHashesScore\n'
                    '  Feature match score: 1.000000\n'
                    '  Feature mismatch score: -1.000000\n'
                    'Overall matches: 0\n'
                    'Significant matches: 0\n'
                    'Query hash count: 0')

        self.assertEqual(expected, result.print_())

    def testPrintWithoutQueryWithNoMatchesDueToNoFinders(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the read and when there are no matches (in this
        case due to the database having no finders).
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [])
        database = Database(params)
        database.addSubject(read)
        findParams = FindParameters(significanceFraction=0.1)
        result = database.find(read, findParams, storeFullAnalysis=True)

        expected = ('Find parameters:\n'
                    '  Significance method: HashFraction\n'
                    '  Significance fraction: 0.100000\n'
                    '  Score method: MinHashesScore\n'
                    '  Feature match score: 1.000000\n'
                    '  Feature mismatch score: -1.000000\n'
                    'Overall matches: 0\n'
                    'Significant matches: 0\n'
                    'Query hash count: 0')
        self.assertEqual(expected, result.print_(printQuery=False))

    def testPrintNoMatchingSubjects(self):
        """
        Check that the print_ method of a result produces the expected result
        when there are no matches.
        """
        query = AARead('query', 'FRRRFRRRFRFRFRFRFRFRFFRRRFRRRFRRRF')
        params = Parameters([AlphaHelix, BetaStrand], [])
        database = Database(params)
        subject = AARead('subject', 'VICVICV')
        database.addSubject(subject)
        result = database.find(query, storeFullAnalysis=True)

        expected = ('Find parameters:\n'
                    '  Significance method: HashFraction\n'
                    '  Significance fraction: 0.250000\n'
                    '  Score method: MinHashesScore\n'
                    '  Feature match score: 1.000000\n'
                    '  Feature mismatch score: -1.000000\n'
                    'Overall matches: 0\n'
                    'Significant matches: 0\n'
                    'Query hash count: 1')
        self.assertEqual(expected, result.print_(printQuery=False))

    def testPrintOneMatchStoredAnalysis(self):
        """
        Check that the print_ method of a result produces the expected result
        when there is a match and we keep the full analysis.
        """
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        params = Parameters([AlphaHelix], [])
        database = Database(params)
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)
        findParams = FindParameters(significanceFraction=0.1,
                                    scoreMethod='MinHashesScore')
        result = database.find(query, findParams, storeFullAnalysis=True)

        expected = ('Find parameters:\n'
                    '  Significance method: HashFraction\n'
                    '  Significance fraction: 0.100000\n'
                    '  Score method: MinHashesScore\n'
                    '  Feature match score: 1.000000\n'
                    '  Feature mismatch score: -1.000000\n'
                    'Overall matches: 1\n'
                    'Significant matches: 1\n'
                    'Query hash count: 1\n'
                    'Matched subjects:\n'
                    '  Subject 1:\n'
                    '    Title: subject\n'
                    '    Best HSP score: 1.0\n'
                    '    Index in database: 0\n'
                    '    Subject hash count: 1\n'
                    '    Subject/query min hash count: 1\n'
                    '    Significance cutoff: 0.100000\n'
                    '    Number of HSPs: 1\n'
                    '      HSP 1 (bin 38): 1 matching hash, score 1.000000\n'
                    '        Landmark AlphaHelix symbol=A offset=0 len=9 '
                    'detail=2\n'
                    '        Trig point AlphaHelix symbol=A offset=25 '
                    'len=13 detail=3')

        self.assertEqual(expected,
                         result.print_(printQuery=False, printFeatures=True))

    def testPrintOneMatchStoredAnalysisPrintHistogram(self):
        """
        Check that the print_ method of a result produces the expected result
        when there is a match and we keep the full analysis and also ask for
        the histogram to be printed.
        """
        params = Parameters([AminoAcidsLm], [])
        database = Database(params)
        subject = AARead('subject', 'CACACAAACACA')
        database.addSubject(subject)
        query = AARead('query', 'CACACA')
        findParams = FindParameters(significanceFraction=0.1,
                                    scoreMethod='MinHashesScore')
        result = database.find(query, findParams, storeFullAnalysis=True)

        self.maxDiff = None
        expected = ('Find parameters:\n'
                    '  Significance method: HashFraction\n'
                    '  Significance fraction: 0.100000\n'
                    '  Score method: MinHashesScore\n'
                    '  Feature match score: 1.000000\n'
                    '  Feature mismatch score: -1.000000\n'
                    'Overall matches: 1\n'
                    'Significant matches: 1\n'
                    'Query hash count: 3\n'
                    'Matched subjects:\n'
                    '  Subject 1:\n'
                    '    Title: subject\n'
                    '    Best HSP score: 1.0\n'
                    '    Index in database: 0\n'
                    '    Subject hash count: 10\n'
                    '    Subject/query min hash count: 3\n'
                    '    Significance cutoff: 0.300000\n'
                    '    Number of HSPs: 6\n'
                    '      HSP 1 (bin 2): 3 matching hashes, score 1.000000\n'
                    '      HSP 2 (bin 0): 1 matching hash, score 0.333333\n'
                    '      HSP 3 (bin 5): 1 matching hash, score 0.333333\n'
                    '      HSP 4 (bin 7): 1 matching hash, score 0.333333\n'
                    '      HSP 5 (bin 10): 1 matching hash, score 0.333333\n'
                    '      HSP 6 (bin 12): 1 matching hash, score 0.333333\n'
                    '    Histogram:\n'
                    '      Number of bins: 13\n'
                    '      Bin width: 0.7692307692\n'
                    '      Max bin count: 3\n'
                    '      Max (scaled) offset delta: 8\n'
                    '      Min (scaled) offset delta: -2\n'
                    '      Non-empty bins:\n'
                    '        Index Count        Range Significant\n'
                    '            0     1 -2.0 to -1.2         Yes\n'
                    '            2     3 -0.5 to +0.3         Yes\n'
                    '            5     1 +1.8 to +2.6         Yes\n'
                    '            7     1 +3.4 to +4.2         Yes\n'
                    '           10     1 +5.7 to +6.5         Yes\n'
                    '           12     1 +7.2 to +8.0         Yes')

        self.assertEqual(expected,
                         result.print_(printQuery=False, printHistograms=True,
                                       printFeatures=False))

    def testOneMatchNoFullAnalysis(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, and when sequences are not printed.
        """
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        params = Parameters([AlphaHelix], [])
        database = Database(params)
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)
        findParams = FindParameters(significanceFraction=0.1,
                                    scoreMethod='MinHashesScore')
        result = database.find(query, findParams)

        expected = ('Find parameters:\n'
                    '  Significance method: HashFraction\n'
                    '  Significance fraction: 0.100000\n'
                    '  Score method: MinHashesScore\n'
                    '  Feature match score: 1.000000\n'
                    '  Feature mismatch score: -1.000000\n'
                    'Overall matches: 1\n'
                    'Significant matches: 1\n'
                    'Query hash count: 1\n'
                    'Matched subjects:\n'
                    '  Subject 1:\n'
                    '    Title: subject\n'
                    '    Best HSP score: 1.0\n'
                    '    Index in database: 0\n'
                    '    Subject hash count: 1\n'
                    '    Subject/query min hash count: 1\n'
                    '    Significance cutoff: 0.100000\n'
                    '    Number of HSPs: 1\n'
                    '      HSP 1 (bin 38): 1 matching hash, score 1.000000\n'
                    '        Landmark AlphaHelix symbol=A offset=0 len=9 '
                    'detail=2\n'
                    '        Trig point AlphaHelix symbol=A offset=25 '
                    'len=13 detail=3')

        self.assertEqual(expected,
                         result.print_(printQuery=False, printFeatures=True))

    def testOneMatchNoFullAnalysisNoFeatures(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, and when features are not printed.
        """
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        params = Parameters([AlphaHelix], [])
        database = Database(params)
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)
        findParams = FindParameters(significanceFraction=0.1,
                                    scoreMethod='MinHashesScore')
        result = database.find(query, findParams)

        expected = ('Find parameters:\n'
                    '  Significance method: HashFraction\n'
                    '  Significance fraction: 0.100000\n'
                    '  Score method: MinHashesScore\n'
                    '  Feature match score: 1.000000\n'
                    '  Feature mismatch score: -1.000000\n'
                    'Overall matches: 1\n'
                    'Significant matches: 1\n'
                    'Query hash count: 1\n'
                    'Matched subjects:\n'
                    '  Subject 1:\n'
                    '    Title: subject\n'
                    '    Best HSP score: 1.0\n'
                    '    Index in database: 0\n'
                    '    Subject hash count: 1\n'
                    '    Subject/query min hash count: 1\n'
                    '    Significance cutoff: 0.100000\n'
                    '    Number of HSPs: 1\n'
                    '      HSP 1 (bin 38): 1 matching hash, score 1.000000')

        self.assertEqual(expected,
                         result.print_(printQuery=False, printFeatures=False))

    def testOneMatchNoFullAnalysisNoFeaturesSixHSPs(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, when features are not printed, and when
        there are six HSPs.
        """
        params = Parameters([AminoAcidsLm], [])
        database = Database(params)
        subject = AARead('subject', 'CACACAAACACA')
        database.addSubject(subject)
        query = AARead('query', 'CACACA')
        findParams = FindParameters(significanceFraction=0.1,
                                    scoreMethod='MinHashesScore')
        result = database.find(query, findParams)

        expected = ('Find parameters:\n'
                    '  Significance method: HashFraction\n'
                    '  Significance fraction: 0.100000\n'
                    '  Score method: MinHashesScore\n'
                    '  Feature match score: 1.000000\n'
                    '  Feature mismatch score: -1.000000\n'
                    'Overall matches: 1\n'
                    'Significant matches: 1\n'
                    'Query hash count: 3\n'
                    'Matched subjects:\n'
                    '  Subject 1:\n'
                    '    Title: subject\n'
                    '    Best HSP score: 1.0\n'
                    '    Index in database: 0\n'
                    '    Subject hash count: 10\n'
                    '    Subject/query min hash count: 3\n'
                    '    Significance cutoff: 0.300000\n'
                    '    Number of HSPs: 6\n'
                    '      HSP 1 (bin 2): 3 matching hashes, score 1.000000\n'
                    '      HSP 2 (bin 0): 1 matching hash, score 0.333333\n'
                    '      HSP 3 (bin 5): 1 matching hash, score 0.333333\n'
                    '      HSP 4 (bin 7): 1 matching hash, score 0.333333\n'
                    '      HSP 5 (bin 10): 1 matching hash, score 0.333333\n'
                    '      HSP 6 (bin 12): 1 matching hash, score 0.333333')

        self.assertEqual(expected,
                         result.print_(printQuery=False, printFeatures=False))

    def testOneMatchNoFullAnalysisNoFeaturesSixHSPsNotSortedByScore(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, when features are not printed, and when
        there are six HSPs and these are not sorted by score (but by histogram
        bin index).
        """
        params = Parameters([AminoAcidsLm], [])
        database = Database(params)
        subject = AARead('subject', 'CACACAAACACA')
        database.addSubject(subject)
        query = AARead('query', 'CACACA')
        findParams = FindParameters(significanceFraction=0.1,
                                    scoreMethod='MinHashesScore')
        result = database.find(query, findParams)

        expected = ('Find parameters:\n'
                    '  Significance method: HashFraction\n'
                    '  Significance fraction: 0.100000\n'
                    '  Score method: MinHashesScore\n'
                    '  Feature match score: 1.000000\n'
                    '  Feature mismatch score: -1.000000\n'
                    'Overall matches: 1\n'
                    'Significant matches: 1\n'
                    'Query hash count: 3\n'
                    'Matched subjects:\n'
                    '  Subject 1:\n'
                    '    Title: subject\n'
                    '    Best HSP score: 1.0\n'
                    '    Index in database: 0\n'
                    '    Subject hash count: 10\n'
                    '    Subject/query min hash count: 3\n'
                    '    Significance cutoff: 0.300000\n'
                    '    Number of HSPs: 6\n'
                    '      HSP 1 (bin 0): 1 matching hash, score 0.333333\n'
                    '      HSP 2 (bin 2): 3 matching hashes, score 1.000000\n'
                    '      HSP 3 (bin 5): 1 matching hash, score 0.333333\n'
                    '      HSP 4 (bin 7): 1 matching hash, score 0.333333\n'
                    '      HSP 5 (bin 10): 1 matching hash, score 0.333333\n'
                    '      HSP 6 (bin 12): 1 matching hash, score 0.333333')

        self.assertEqual(expected,
                         result.print_(printQuery=False, printFeatures=False,
                                       sortHSPsByScore=False))

    def testPrintWithoutQueryWithMatchesWithoutFullAnalysisWithSequences(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, and when sequences are printed.
        """
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        params = Parameters([AlphaHelix], [])
        database = Database(params)
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)
        findParams = FindParameters(significanceFraction=0.1,
                                    scoreMethod='MinHashesScore')
        result = database.find(query, findParams)

        expected = (
            'Find parameters:\n'
            '  Significance method: HashFraction\n'
            '  Significance fraction: 0.100000\n'
            '  Score method: MinHashesScore\n'
            '  Feature match score: 1.000000\n'
            '  Feature mismatch score: -1.000000\n'
            'Overall matches: 1\n'
            'Significant matches: 1\n'
            'Query hash count: 1\n'
            'Matched subjects:\n'
            '  Subject 1:\n'
            '    Title: subject\n'
            '    Best HSP score: 1.0\n'
            '    Sequence: FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF\n'
            '    Index in database: 0\n'
            '    Subject hash count: 1\n'
            '    Subject/query min hash count: 1\n'
            '    Significance cutoff: 0.100000\n'
            '    Number of HSPs: 1\n'
            '      HSP 1 (bin 38): 1 matching hash, score '
            '1.000000')

        self.assertEqual(expected,
                         result.print_(printQuery=False, printSequences=True,
                                       printFeatures=False))

    def testBinOverflowWarning(self):
        """
        When a large distanceBase is used, it scales many distances into the
        same value. This can cause more deltas to be placed into the same
        histogram bin than there are hashes in the query. Normally, the
        maximum number of deltas in a bin will be the number of hashes in
        the query, so anything more than this indicates that the
        distanceBase is (in this case) clouding the picture because it is
        too aggressive. That's not necessarily a bad thing - a high
        distanceBase might be good overall in some circumstances (e.g.,
        comparing homologous sequences that have had many insertions /
        deletions), and just undesirable in one particular bin. So for now
        a warning is raised in result.py. In this test we trigger the
        warning and catch it.

        When this condition occurs, the score from find() must be set to 1.0
        """
        # Here's a sequence with an alpha helix in the middle of many troughs.
        sequence = (
            'ISDLFKKNQHGGLREIYVLTIKSKLLALFLETCSRCLCEQFTVETMTHPDCKMEVIERHKMQVN'
            'TLAAKTKRSFNSYQCSADKKSWNNNLVMPALSIPLLMLLPKAFHGAIQRTLNMWNQRLIQLPRG'
            'FRRRFRRRF'
            'VLKLLTAGVELSDPTYQLMMKEFESPGSTGDSPLFPNAGSGYCILHRGMMQGILHYTSSLLHVN'
            'YLFVTRELIRSAYKAKFPDTTFLIDQMCSSDDSATIMSVVHPLNESEQGIKVISAFSEIICEVL'
            'KTFCRYSCFTNSEKSVMGSLNQLEFNSEFIIGNNMAVPILKWVFSAFG')
        params = Parameters([AlphaHelix], [Troughs], distanceBase=1.4)
        database = Database(params)
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            # This find will trigger a warning due to too many deltas.
            result = database.find(query)
            self.assertEqual(1, len(w))
            self.assertTrue(issubclass(w[0].category, RuntimeWarning))
            error = ('Bin contains 14 deltas for a query/subject pair with a '
                     'minimum hash count of only 10.')
            self.assertIn(error, str(w[0].message))
            self.assertEqual(1.0, result.analysis['0']['bestScore'])

    def testCorrectSignificanceAnalysisAlways(self):
        """
        When storeFullAnalysis is True, the result must have the correct
        significance analysis if the Always significance method was used.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([], [])
        database = Database(params)
        database.addSubject(AARead('subject', 'A' * 21))
        hashCount = 1
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 1, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 2, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='Always',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.1)
        result = Result(read, database, matches, hashCount, findParams,
                        storeFullAnalysis=True)
        significanceAnalysis = result.analysis['0']['significanceAnalysis']
        self.assertEqual('Always', significanceAnalysis['significanceMethod'])

    def testScoreAnalysisAddedToResult(self):
        """
        When storeFullAnalysis is True, the result must have a scoreAnalysis
        in each significant bin.
        """
        read = AARead('read', 'AGTARFSDDD')
        params = Parameters([AlphaHelix], [AminoAcids])
        database = Database(params)
        _, subjectIndex1, _ = database.addSubject(AARead('subject0',
                                                         'ADDDADDDAWWWWW'))
        _, subjectIndex2, _ = database.addSubject(AARead('subject1',
                                                         'ADDDADDDAWWWWW'))
        hashCount = 5
        matches = {
            '0': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 2),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 3),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 10, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],

            '1': [
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 1),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
                {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 2),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 0),
                },
            ],
        }
        findParams = FindParameters(significanceMethod='HashFraction',
                                    scoreMethod='MinHashesScore',
                                    significanceFraction=0.3)
        result = Result(read, database, matches, hashCount, findParams,
                        storeFullAnalysis=True)

        # Note that the contents of the score analysis for each scoring
        # method are checked in test/test_score.py, so here we just check
        # that the Result class has made the analysis available.
        for subjectIndex in subjectIndex1, subjectIndex2:
            self.assertIn('scoreAnalysis',
                          result.analysis[subjectIndex]['significantBins'][0])
