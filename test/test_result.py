from unittest import TestCase
from cStringIO import StringIO
from json import loads
import warnings

from dark.reads import AARead

from light.result import Result
from light.features import Landmark, TrigPoint
from light.reads import ScannedRead
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
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        hashCount = 0
        result = Result(read, {}, hashCount, significanceFraction=0.1,
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
        result = Result(read, matches, hashCount, significanceFraction=0.1,
                        database=database)
        self.assertEqual(matches, result.matches)

    def testNoSignificantMatches(self):
        """
        No matches are significant if there are not enough distance deltas
        above the mean number of deltas.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([AlphaHelix], [AminoAcids])
        database.addSubject(AARead('subject', 'ADDDADDDAWW'))
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
        database = Database([AlphaHelix], [AminoAcids])
        database.addSubject(AARead('subject', 'ADDDADDDAWWWW'))
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
        database = Database([AlphaHelix], [AminoAcids])
        # Note that both subject added here have a hash count of 5 (the
        # same as the passed query hash count). If we use anything less,
        # the score will change because its calculation uses the min
        # query/subject hash count in the denominator.
        database.addSubject(AARead('subject0', 'ADDDADDDAWWWWW'))
        database.addSubject(AARead('subject1', 'ADDDADDDAWWWWW'))
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
        result = Result(read, {}, 0, significanceFraction=0.1,
                        database=database)
        fp = StringIO()
        self.assertIs(fp, result.save(fp))

    def testSave(self):
        """
        Save must write results out in the expected JSON format.
        """
        read = ScannedRead(
            AARead('id', 'FRRRFRRRFRFRFRFRRRFRRRF'))
        database = Database([AlphaHelix], [])
        database.addSubject(AARead('subject0', 'FRRRFRRRFRFRFRFRRRFRRRF'))
        database.addSubject(AARead('subject1', 'FRRRFRRRFRFRFRFRRRFRRRF'))
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
        self.maxDiff = None
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
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        database.addSubject(AARead('subject', 'A' * 21))
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
        result = Result(read, matches, hashCount, significanceFraction=0.1,
                        database=database, storeFullAnalysis=True)
        self.assertEqual(21, result.analysis[0]['histogram'].nBins)

    def testRightNumberOfBucketsDefaultNonEven(self):
        """
        If no distanceBase is specified for a database, the number of bins must
        be 21 if the length of the subject is 20 and the length of the read
        is less. This is because int(log base 1.1 of 20) = 31, which is greater
        than 20, so the value of 20 should be used but result.py will adjust
        this to 21 to ensure that the number of bins is odd.
        """
        read = ScannedRead(AARead('read', 'AGTARFSDDD'))
        database = Database([], [])
        database.addSubject(AARead('subject', 'A' * 20))
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
        result = Result(read, matches, hashCount, significanceFraction=0.1,
                        database=database, storeFullAnalysis=True)
        self.assertEqual(21, result.analysis[0]['histogram'].nBins)

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
        result = Result(read, matches, hashCount, significanceFraction=0.1,
                        database=database, storeFullAnalysis=True)
        self.assertEqual(11, result.analysis[0]['histogram'].nBins)

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
        result = database.find(read, significanceFraction=0.1,
                               storeFullAnalysis=True)

        result.print_(fp=fp)
        expected = ("Query title: read\n"
                    "  Length: 10\n"
                    "  Covered indices: 0 (0.00%)\n"
                    "  Landmark count 0, trig point count 0\n"
                    "Overall matches: 0\n"
                    "Significant matches: 0\n"
                    "Query hash count: 0\n"
                    "Significance fraction: 0.100000\n")

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
        result = database.find(read, significanceFraction=0.1,
                               storeFullAnalysis=True)

        result.print_(fp=fp, printQuery=False)
        expected = ("Overall matches: 0\n"
                    "Significant matches: 0\n"
                    "Query hash count: 0\n"
                    "Significance fraction: 0.100000\n")
        self.assertEqual(expected, fp.getvalue())

    def testPrintNoMatchingSubjects(self):
        """
        Check that the print_ method of a result produces the expected result
        when there are no matches.
        """
        fp = StringIO()
        query = AARead('query', 'FRRRFRRRFRFRFRFRFRFRFFRRRFRRRFRRRF')
        database = Database([AlphaHelix, BetaStrand], [])
        subject = AARead('subject', 'VICVICV')
        database.addSubject(subject)
        result = database.find(query, storeFullAnalysis=True)

        result.print_(fp=fp, printQuery=False)
        expected = ("Overall matches: 0\n"
                    "Significant matches: 0\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.250000\n")
        self.assertEqual(expected, fp.getvalue())

    def testPrintOneMatchStoredAnalysis(self):
        """
        Check that the print_ method of a result produces the expected result
        when there is a match and we keep the full analysis.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)
        result = database.find(query, significanceFraction=0.1,
                               storeFullAnalysis=True)

        result.print_(fp=fp, printQuery=False)

        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.100000\n"
                    "Matched subjects:\n"
                    "  Subject 1:\n"
                    "    Title: subject\n"
                    "    Best HSP score: 1.0\n"
                    "    Index in database: 0\n"
                    "    Subject hash count: 1\n"
                    "    Subject/query min hash count: 1\n"
                    "    Significance cutoff: 0.100000\n"
                    "    Number of HSPs: 1\n"
                    "      HSP 1 (bin 38): 1 matching hash, score 1.000000\n"
                    "        Landmark AlphaHelix symbol='A' offset=0 len=9 "
                    "detail=2 subjectOffset=0\n"
                    "        Trig point AlphaHelix symbol='A' offset=25 "
                    "len=13 detail=3\n")

        self.assertEqual(expected, fp.getvalue())

    def testPrintOneMatchStoredAnalysisPrintHistogram(self):
        """
        Check that the print_ method of a result produces the expected result
        when there is a match and we keep the full analysis and also ask for
        the histogram to be printed.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)
        result = database.find(query, significanceFraction=0.1,
                               storeFullAnalysis=True)

        result.print_(fp=fp, printQuery=False, printHistograms=True)

        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.100000\n"
                    "Matched subjects:\n"
                    "  Subject 1:\n"
                    "    Title: subject\n"
                    "    Best HSP score: 1.0\n"
                    "    Index in database: 0\n"
                    "    Subject hash count: 1\n"
                    "    Subject/query min hash count: 1\n"
                    "    Significance cutoff: 0.100000\n"
                    "    Number of HSPs: 1\n"
                    "      HSP 1 (bin 38): 1 matching hash, score 1.000000\n"
                    "        Landmark AlphaHelix symbol='A' offset=0 len=9 "
                    "detail=2 subjectOffset=0\n"
                    "        Trig point AlphaHelix symbol='A' offset=25 "
                    "len=13 detail=3\n"
                    "    Histogram:\n"
                    "      Number of bins: 39\n"
                    "      Bin width: 0.0000000000\n"
                    "      Max bin count: 1\n"
                    "      Non-empty bins (index:count; *=significant): "
                    "38:1*\n"
                    "      Max (scaled) offset delta: 0\n"
                    "      Min (scaled) offset delta: 0\n")

        self.assertEqual(expected, fp.getvalue())

    def testOneMatchNoFullAnalysis(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, and when sequences are not printed.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)
        result = database.find(query, significanceFraction=0.1)

        result.print_(fp=fp, printQuery=False)
        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.100000\n"
                    "Matched subjects:\n"
                    "  Subject 1:\n"
                    "    Title: subject\n"
                    "    Best HSP score: 1.0\n"
                    "    Index in database: 0\n"
                    "    Subject hash count: 1\n"
                    "    Subject/query min hash count: 1\n"
                    "    Significance cutoff: 0.100000\n"
                    "    Number of HSPs: 1\n"
                    "      HSP 1 (bin 38): 1 matching hash, score 1.000000\n"
                    "        Landmark AlphaHelix symbol='A' offset=0 len=9 "
                    "detail=2 subjectOffset=0\n"
                    "        Trig point AlphaHelix symbol='A' offset=25 "
                    "len=13 detail=3\n")

        self.assertEqual(expected, fp.getvalue())

    def testOneMatchNoFullAnalysisNoFeatures(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, and when features are not printed.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)
        result = database.find(query, significanceFraction=0.1)

        result.print_(fp=fp, printQuery=False, printFeatures=False)
        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.100000\n"
                    "Matched subjects:\n"
                    "  Subject 1:\n"
                    "    Title: subject\n"
                    "    Best HSP score: 1.0\n"
                    "    Index in database: 0\n"
                    "    Subject hash count: 1\n"
                    "    Subject/query min hash count: 1\n"
                    "    Significance cutoff: 0.100000\n"
                    "    Number of HSPs: 1\n"
                    "      HSP 1 (bin 38): 1 matching hash, score 1.000000\n")

        self.assertEqual(expected, fp.getvalue())

    def testOneMatchNoFullAnalysisNoFeaturesSixHSPs(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, when features are not printed, and when
        there are six HSPs.
        """
        fp = StringIO()
        subject = 'CACACAAACACA'
        query = 'CACACA'
        database = Database([AminoAcidsLm], [])
        subject = AARead('subject', subject)
        database.addSubject(subject)
        query = AARead('query', query)
        result = database.find(query, significanceFraction=0.1)

        result.print_(fp=fp, printQuery=False, printFeatures=False)
        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 3\n"
                    "Significance fraction: 0.100000\n"
                    "Matched subjects:\n"
                    "  Subject 1:\n"
                    "    Title: subject\n"
                    "    Best HSP score: 1.0\n"
                    "    Index in database: 0\n"
                    "    Subject hash count: 10\n"
                    "    Subject/query min hash count: 3\n"
                    "    Significance cutoff: 0.300000\n"
                    "    Number of HSPs: 6\n"
                    "      HSP 1 (bin 2): 3 matching hashes, score 1.000000\n"
                    "      HSP 2 (bin 0): 1 matching hash, score 0.333333\n"
                    "      HSP 3 (bin 5): 1 matching hash, score 0.333333\n"
                    "      HSP 4 (bin 7): 1 matching hash, score 0.333333\n"
                    "      HSP 5 (bin 10): 1 matching hash, score 0.333333\n"
                    "      HSP 6 (bin 12): 1 matching hash, score 0.333333\n")

        self.assertEqual(expected, fp.getvalue())

    def testOneMatchNoFullAnalysisNoFeaturesSixHSPsNotSortedByScore(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, when features are not printed, and when
        there are six HSPs and these are not sorted by score (but by histogram
        bin index).
        """
        fp = StringIO()
        subject = 'CACACAAACACA'
        query = 'CACACA'
        database = Database([AminoAcidsLm], [])
        subject = AARead('subject', subject)
        database.addSubject(subject)
        query = AARead('query', query)
        result = database.find(query, significanceFraction=0.1)

        result.print_(fp=fp, printQuery=False, printFeatures=False,
                      sortHSPsByScore=False)
        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 3\n"
                    "Significance fraction: 0.100000\n"
                    "Matched subjects:\n"
                    "  Subject 1:\n"
                    "    Title: subject\n"
                    "    Best HSP score: 1.0\n"
                    "    Index in database: 0\n"
                    "    Subject hash count: 10\n"
                    "    Subject/query min hash count: 3\n"
                    "    Significance cutoff: 0.300000\n"
                    "    Number of HSPs: 6\n"
                    "      HSP 1 (bin 0): 1 matching hash, score 0.333333\n"
                    "      HSP 2 (bin 2): 3 matching hashes, score 1.000000\n"
                    "      HSP 3 (bin 5): 1 matching hash, score 0.333333\n"
                    "      HSP 4 (bin 7): 1 matching hash, score 0.333333\n"
                    "      HSP 5 (bin 10): 1 matching hash, score 0.333333\n"
                    "      HSP 6 (bin 12): 1 matching hash, score 0.333333\n")

        self.assertEqual(expected, fp.getvalue())

    def testPrintWithoutQueryWithMatchesWithoutFullAnalysisWithSequences(self):
        """
        Check that the print_ method of a result produces the expected result
        when asked to not print the query, when there are matches and the
        full analysis is not stored, and when sequences are printed.
        """
        fp = StringIO()
        sequence = 'FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF'
        database = Database([AlphaHelix], [])
        subject = AARead('subject', sequence)
        database.addSubject(subject)
        query = AARead('query', sequence)
        result = database.find(query, significanceFraction=0.1)

        result.print_(fp=fp, printQuery=False, printSequences=True,
                      printFeatures=False)
        expected = ("Overall matches: 1\n"
                    "Significant matches: 1\n"
                    "Query hash count: 1\n"
                    "Significance fraction: 0.100000\n"
                    "Matched subjects:\n"
                    "  Subject 1:\n"
                    "    Title: subject\n"
                    "    Best HSP score: 1.0\n"
                    "    Sequence: FRRRFRRRFRFRFRFRFRFRFRFRFFRRRFRRRFRRRF\n"
                    "    Index in database: 0\n"
                    "    Subject hash count: 1\n"
                    "    Subject/query min hash count: 1\n"
                    "    Significance cutoff: 0.100000\n"
                    "    Number of HSPs: 1\n"
                    "      HSP 1 (bin 38): 1 matching hash, score "
                    "1.000000\n")

        self.assertEqual(expected, fp.getvalue())

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
        database = Database([AlphaHelix], [Troughs], distanceBase=1.4)
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
            self.assertEqual(1.0, result.analysis[0]['bestScore'])
