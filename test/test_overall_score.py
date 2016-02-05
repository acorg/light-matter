import six
from unittest import TestCase

from dark.reads import AARead, AAReadWithX

from light.database import DatabaseSpecifier
from light.features import Landmark, TrigPoint
from light.histogram import Histogram
from light.landmarks import AlphaHelix, AminoAcids as AminoAcidsLm
from light.trig.amino_acids import AminoAcids
from light.trig.peaks import Peaks
from light.overall_score import (
    BestBinScore, SignificantBinScore, offsetsInBin, computeLengthNormalizer,
    GreedySignificantBinScore)
from light.parameters import Parameters, FindParameters
from light.subject import Subject


class TestBestBinScore(TestCase):
    """
    Tests for the light.overall_score.BestBinScore class.
    """
    def testEmptyHistogram(self):
        """
        An empty histogram must return a score of C{None}.
        """
        histogram = Histogram()
        histogram.finalize()
        bestBinScore = BestBinScore(histogram, [])
        score, analysis = bestBinScore.calculateScore()
        self.assertEqual(None, score)
        self.assertEqual(
            {
                'score': score,
                'scoreClass': bestBinScore.__class__,
            },
            analysis)

    def testCompareEqualSequencesScoreMustBeOne(self):
        """
        If a sequence is compared to itself, the overall score must be 1.0. See
        https://github.com/acorg/light-matter/issues/321.
        This is a real-life test that it actually works.
        """
        pichninde = AARead('pichninde', 'RLKFGLSYKEQVGGNRELYVGDLNTKLTTRLIEDYS'
                                        'ESLMQNMRYTCLNNEKEFERALLDMKSVVRQSGLAV'
                                        'SMDHSKWGPHMSPVIFAALLKGLEFNLKDGSEVPNA'
                                        'AIVNILLWHIHKMVEVPFNVVEAYMKGFLKRGLGMM'
                                        'DKGGCTIAEEFMFGYFEKGKVPSHISSVLDMGQGIL'
                                        'HNTSDLYGLITEQFINYALELCYGVRFISYTSSDDE'
                                        'IMLSLNEAFKFKDRDELNVDLVLDCMEFHYFLSDKL'
                                        'NKFVSPKTVVGTFASEFKSRFFIWSQEVPLLTKFVA'
                                        'AALH')

        db = DatabaseSpecifier().getDatabaseFromKeywords(
            landmarkNames=[
                'AlphaHelix', 'AlphaHelix_3_10', 'AlphaHelix_pi',
                'AminoAcidsLm', 'BetaStrand', 'BetaTurn', 'Prosite'],
            trigPointNames=['AminoAcids', 'Peaks', 'Troughs'],
            distanceBase=1.01, limitPerLandmark=50, minDistance=1,
            maxDistance=100)
        _, subjectIndex, _ = db.addSubject(pichninde)

        findParams = FindParameters(significanceFraction=0.01,
                                    scoreMethod='FeatureAAScore',
                                    overallScoreMethod='BestBinScore')
        result = db.find(pichninde, findParams, storeFullAnalysis=True)
        self.assertEqual(1.0, result.analysis[subjectIndex]['overallScore'])

    def testCompareRealSequencesBestBinScoreAndOverallScoreIdentical(self):
        """
        When comparing two sequences, the best bin score (calculated using the
        FeatureAAScore class) and the overall score (calculated using the
        BestBinScore class which sets the overall score to the best bin score)
        must be the same.
        """
        golv = AARead('GOLV', 'RVDIFKKNQHGGLREIYVLDLASRIVQLCLEEISRAVCQELPIEMM'
                              'MHPELKLKKPQEHMYKAAISPESYKSNVSSSNDAKVWNQGHHVAKF'
                              'AQFLCRLLSPEWHGLIVNGLKLWTNKKIALPDGVMNILSRANTPLF'
                              'RNSIHQAVHDSYKGITPMRWLRPGETFMRIESGMMQGILHYTSSLF'
                              'HASLLMMRDSLWRSYSEQLGVKSITTDLVSSDDSSRMTDIFYRDSK'
                              'NFKRGKIFARADHMAIEPLSRCFGIWMSPKSTYCCNGIMEFNSEYF'
                              'FRASLYRPTLKWSYACLG')

        akav = AARead('AKAV', 'VFTYFNKGQKTAKDREIFVGEFEAKMCLYLVERISKERCKLNPDEM'
                              'ISEPGDGKLKKLEDMAEYEIRYTANTLKSMKDKALQEFSKFADDFN'
                              'FKPHSTKIEINADMSKWSAQDVLFKYFWLFALDPALYKPEKERILY'
                              'FLCNYMDKVLVIPDDVMTSILDQRVKREKDIIYEMTNGLKQNWVSI'
                              'KRNWLQGNLNYTSSYLHSCCMNVYKDIIKNVATLLEGDVLVNSMVH'
                              'SDDNHTSITMIQDKFPDDIIIEYCIKLFEKICLSFGNQANMKKTYV'
                              'TNFIKEFVSLFNIYGEPFSVYGRFLLTAVG')

        findParams = FindParameters(significanceFraction=0.01,
                                    scoreMethod='FeatureAAScore')

        kwds = dict(landmarkNames=['Prosite'], trigPointNames=['Peaks'],
                    distanceBase=1, limitPerLandmark=40, minDistance=1,
                    maxDistance=10000)

        db = DatabaseSpecifier().getDatabaseFromKeywords(**kwds)
        _, subjectIndex, _ = db.addSubject(golv)
        result = db.find(akav, findParams, storeFullAnalysis=True)
        self.assertEqual(0.27265479670475684,
                         result.analysis[subjectIndex]['overallScore'])
        self.assertEqual(result.analysis[subjectIndex]['bestBinScore'],
                         result.analysis[subjectIndex]['overallScore'])

    def testPrintAnalysis(self):
        """
        The analysis of the overall score calculation must print correctly.
        """
        histogram = Histogram()
        histogram.finalize()
        bestBinScore = BestBinScore(histogram, [])
        score, analysis = bestBinScore.calculateScore()

        expected = (
            'Overall score method: BestBinScore\nOverall score: %s' % score)

        self.assertEqual(expected, BestBinScore.printAnalysis(analysis))


class TestOffsetsInBin(TestCase):
    """
    Tests for the light.overall_score.offsetsInBin function.
    """
    def testIncorrectQueryOrSubject(self):
        """
        If the value of queryOrSubject passed to offsetsInBin is not
        'query' or 'subject', a C{KeyError} must be raised.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        error = 'xxxLandmark'
        six.assertRaisesRegex(self, KeyError, error, offsetsInBin,
                              histogram.bins[0], 'xxx',
                              {queryLandmark, queryTrigPoint})

    def testQueryOffsetsOnePairInBin(self):
        """
        The query offets in a bin must be found correctly when there is one
        pair in the bin and no other features exist.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        matchedOffsets, unmatchedOffsets, minOffset, maxOffset = (
            offsetsInBin(histogram.bins[0], 'query',
                         {queryLandmark, queryTrigPoint}))
        self.assertEqual(set(range(100, 120)), matchedOffsets)
        self.assertEqual(set(), unmatchedOffsets)
        self.assertEqual(minOffset, 100)
        self.assertEqual(maxOffset, 119)

    def testSubjectOffsetsOnePairInBin(self):
        """
        The subject offets in a bin must be found correctly when there is one
        pair in the bin and no other features exist.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 200, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 210)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        matchedOffsets, unmatchedOffsets, minOffset, maxOffset = (
            offsetsInBin(histogram.bins[0], 'subject',
                         {subjectLandmark, subjectTrigPoint}))
        self.assertEqual(set(range(200, 220)), matchedOffsets)
        self.assertEqual(set(), unmatchedOffsets)
        self.assertEqual(minOffset, 200)
        self.assertEqual(maxOffset, 219)

    def testQueryOffsetsTwoPairsInBin(self):
        """
        The query offets in a bin must be found correctly when there are two
        pairs in the bin and no other features exist.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 110)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 210)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 210)

        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark1,
            'queryTrigPoint': queryTrigPoint1,
            'subjectLandmark': subjectLandmark1,
            'subjectTrigPoint': subjectTrigPoint1,
        })
        histogram.add(44, {
            'queryLandmark': queryLandmark2,
            'queryTrigPoint': queryTrigPoint2,
            'subjectLandmark': subjectLandmark2,
            'subjectTrigPoint': subjectTrigPoint2,
        })
        histogram.finalize()
        matchedOffsets, unmatchedOffsets, minOffset, maxOffset = (
            offsetsInBin(histogram.bins[0], 'query',
                         {queryLandmark1, queryTrigPoint1,
                          queryLandmark2, queryTrigPoint2}))
        self.assertEqual(set(range(100, 120)) | set(range(200, 220)),
                         matchedOffsets)
        self.assertEqual(set(), unmatchedOffsets)
        self.assertEqual(minOffset, 100)
        self.assertEqual(maxOffset, 219)

    def testQueryOffsetsTwoPairsInBinUnmatchedFeatureBetween(self):
        """
        The query offets in a bin must be found correctly when there are two
        pairs in the bin and one unmatched feature exists between the pairs.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 110)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 210)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 210)

        unmatchedLandmark = Landmark('AlphaHelix', 'A', 150, 10)
        unmatchedTrigPoint = TrigPoint('Peaks', 'P', 170)

        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark1,
            'queryTrigPoint': queryTrigPoint1,
            'subjectLandmark': subjectLandmark1,
            'subjectTrigPoint': subjectTrigPoint1,
        })
        histogram.add(44, {
            'queryLandmark': queryLandmark2,
            'queryTrigPoint': queryTrigPoint2,
            'subjectLandmark': subjectLandmark2,
            'subjectTrigPoint': subjectTrigPoint2,
        })
        histogram.finalize()
        matchedOffsets, unmatchedOffsets, minOffset, maxOffset = (
            offsetsInBin(histogram.bins[0], 'query',
                         {queryLandmark1, queryTrigPoint1,
                          queryLandmark2, queryTrigPoint2,
                          unmatchedLandmark, unmatchedTrigPoint}))
        self.assertEqual(set(range(100, 120)) | set(range(200, 220)),
                         matchedOffsets)
        self.assertEqual(set(range(150, 160)) | {170}, unmatchedOffsets)
        self.assertEqual(minOffset, 100)
        self.assertEqual(maxOffset, 219)

    def testOnePairInBinOneUnmatchedFeatureOverlapped(self):
        """
        The query offets in a bin must be found correctly when there is one
        pair in the bin and one other unmatched feature, but the unmatched
        feature is overlapped by a matching feature. In this case the unmatched
        offsets is empty.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 110)
        unmatchedLandmark = Landmark('AlphaHelix', 'A', 105, 10)
        unmatchedTrigPoint = TrigPoint('Peaks', 'P', 118)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        matchedOffsets, unmatchedOffsets, minOffset, maxOffset = (
            offsetsInBin(histogram.bins[0], 'query',
                         {queryLandmark, queryTrigPoint,
                          unmatchedLandmark, unmatchedTrigPoint}))
        self.assertEqual(set(range(100, 120)), matchedOffsets)
        self.assertEqual(set(), unmatchedOffsets)
        self.assertEqual(minOffset, 100)
        self.assertEqual(maxOffset, 119)

    def testOnePairInBinOneUnmatchedFeatureNotOverlapped(self):
        """
        The query offets in a bin must be found correctly when there is one
        pair in the bin and one other unmatched feature, and the unmatched
        feature is not overlapped by a matching feature. In this case the
        unmatched offsets will contain the offsets of the unmatched feature.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 150)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 150)
        unmatchedLandmark = Landmark('AlphaHelix', 'A', 130, 10)
        unmatchedTrigPoint = TrigPoint('Peaks', 'P', 145)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        matchedOffsets, unmatchedOffsets, minOffset, maxOffset = (
            offsetsInBin(histogram.bins[0], 'query',
                         {queryLandmark, queryTrigPoint,
                          unmatchedLandmark, unmatchedTrigPoint}))
        self.assertEqual(set(range(100, 120)) | {150}, matchedOffsets)
        self.assertEqual(set(range(130, 140)) | {145}, unmatchedOffsets)
        self.assertEqual(minOffset, 100)
        self.assertEqual(maxOffset, 150)

    def testOnePairInBinOneUnmatchedFeaturePartlyOverlapped(self):
        """
        The query offets in a bin must be found correctly when there is one
        pair in the bin and one other unmatched feature, and the unmatched
        feature is partially overlapped by a matching feature. In this case the
        unmatched offsets will contain the non-overlapped offsets of the
        unmatched feature.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 150)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 150)
        unmatchedLandmark = Landmark('AlphaHelix', 'A', 115, 10)
        unmatchedTrigPoint = TrigPoint('Peaks', 'P', 145)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        matchedOffsets, unmatchedOffsets, minOffset, maxOffset = (
            offsetsInBin(histogram.bins[0], 'query',
                         {queryLandmark, queryTrigPoint,
                          unmatchedLandmark, unmatchedTrigPoint}))
        self.assertEqual(set(range(100, 120)) | {150}, matchedOffsets)
        self.assertEqual(set(range(120, 125)) | {145}, unmatchedOffsets)
        self.assertEqual(minOffset, 100)
        self.assertEqual(maxOffset, 150)


class TestComputeLengthNormalizer(TestCase):
    """
    Tests for the light.overall_score.computeLengthNormalizer function.
    """
    def testEmptyArguments(self):
        """
        If there are no features or offsets, the normalizer must be 1.0
        and the numerator and denominator both 0.0.
        """
        normalizer, numerator, denominator = computeLengthNormalizer(
            {}, {}, {}, {})
        self.assertEqual(normalizer, 1.0)
        self.assertEqual(numerator, 0)
        self.assertEqual(denominator, 0)

    def testNoFeaturesMatchedRegion(self):
        """
        If there are no features outside the matched region, the normalizer
        must be 1.0 and the numerator and denominator must both contain the
        number of offsets in features in the matched region (whether or not
        the features are part of the match).
        """
        matchedLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        matchedTrigPoint = TrigPoint('Peaks', 'P', 150)
        unmatchedLandmark = Landmark('AlphaHelix', 'A', 200, 10)
        unmatchedTrigPoint = TrigPoint('Peaks', 'P', 250)

        allFeatures = {matchedLandmark, matchedTrigPoint,
                       unmatchedLandmark, unmatchedTrigPoint}
        overallMatchedOffsets = set(range(100, 120)) | {150}
        overallUnmatchedOffsets = set(range(200, 210)) | {250}
        offsetsInBins = overallMatchedOffsets | overallUnmatchedOffsets

        normalizer, numerator, denominator = computeLengthNormalizer(
            allFeatures, overallMatchedOffsets, overallUnmatchedOffsets,
            offsetsInBins)

        self.assertEqual(normalizer, 1.0)
        self.assertEqual(numerator, 32)
        self.assertEqual(denominator, 32)

    def testOneFeatureOutsideMatchedRegion(self):
        """
        If there is one feature outside the matched region, the normalizer,
        numerator, and denominator must all be as expected.
        """
        matchedLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        matchedTrigPoint = TrigPoint('Peaks', 'P', 150)
        unmatchedLandmark = Landmark('AlphaHelix', 'A', 200, 10)
        unmatchedTrigPoint = TrigPoint('Peaks', 'P', 250)
        otherLandmark = Landmark('AlphaHelix', 'A', 300, 20)
        otherTrigPoint = TrigPoint('Peaks', 'P', 350)

        allFeatures = {matchedLandmark, matchedTrigPoint,
                       unmatchedLandmark, unmatchedTrigPoint,
                       otherLandmark, otherTrigPoint}
        overallMatchedOffsets = set(range(100, 120)) | {150}
        overallUnmatchedOffsets = set(range(200, 210)) | {250}
        offsetsInBins = overallMatchedOffsets | overallUnmatchedOffsets

        normalizer, numerator, denominator = computeLengthNormalizer(
            allFeatures, overallMatchedOffsets, overallUnmatchedOffsets,
            offsetsInBins)

        self.assertEqual(normalizer, 32 / 53)
        self.assertEqual(numerator, 32)
        self.assertEqual(denominator, 53)

    def testOneFeaturePartiallyOutsideMatchedRegion(self):
        """
        If there is one feature partially outside the matched region, the
        normalizer, numerator, and denominator must all be as expected.
        """
        matchedLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        matchedTrigPoint = TrigPoint('Peaks', 'P', 150)
        unmatchedLandmark = Landmark('AlphaHelix', 'A', 200, 10)
        unmatchedTrigPoint = TrigPoint('Peaks', 'P', 250)
        # The first 20 offsets (80-99) of otherLandmark fall outside
        # the matched region.
        otherLandmark = Landmark('AlphaHelix', 'A', 80, 40)
        otherTrigPoint = TrigPoint('Peaks', 'P', 110)

        allFeatures = {matchedLandmark, matchedTrigPoint,
                       unmatchedLandmark, unmatchedTrigPoint,
                       otherLandmark, otherTrigPoint}
        overallMatchedOffsets = set(range(100, 120)) | {150}
        overallUnmatchedOffsets = set(range(200, 210)) | {250}
        offsetsInBins = overallMatchedOffsets | overallUnmatchedOffsets

        normalizer, numerator, denominator = computeLengthNormalizer(
            allFeatures, overallMatchedOffsets, overallUnmatchedOffsets,
            offsetsInBins)

        self.assertEqual(normalizer, 32 / 52)
        self.assertEqual(numerator, 32)
        self.assertEqual(denominator, 52)


class TestSignificantBinScore(TestCase):
    """
    Tests for the light.overall_score.SignificantBinScore class.
    """
    def testNoSignificantBins(self):
        """
        If there are no significant bins, the overall score must be C{None}.
        """
        histogram = Histogram()
        histogram.finalize()
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        sbs = SignificantBinScore([], query, subject, params)
        score, analysis = sbs.calculateScore()
        self.assertIs(None, score)
        self.assertEqual({'score': None,
                          'scoreClass': SignificantBinScore}, analysis)

    def testQueryHasOneFeature(self):
        """
        If the query has one feature, but no pair, the match must have an
        overall score of C{None}.
        """
        histogram = Histogram(1)
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        sbs = SignificantBinScore([], query, subject, params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(None, score)
        self.assertEqual({'score': None,
                          'scoreClass': SignificantBinScore}, analysis)

    def testQueryAndSubjectHaveOneFeature(self):
        """
        If the query and the subject have one feature each, but no pairs, the
        match must have an overall score of C{None}.
        """
        histogram = Histogram(1)
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'AAAAAAAAAAAAAAFRRRFRRRF', 0)
        sbs = SignificantBinScore([], query, subject, params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(None, score)
        self.assertEqual({'score': None,
                          'scoreClass': SignificantBinScore}, analysis)

    def testOnePairInOneBin(self):
        """
        A match with one bin containing one pair must have an overall score of
        1.0 if the query and subject have no additional (non-matching) hashes.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 20,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'overallScoreAdjustedToBestBinScore': False,
                'queryOffsetsInBins': 20,
                'normalizerQuery': 1.0,
                'normalizerSubject': 1.0,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 20,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOnePairInBinOccurringInTwoPlaces(self):
        """
        A match with one bin containing two pairs that occur in two places must
        have an overall score of 1.0 if the query and subject have no
        additional (non-matching) hashes.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 0, 9)
        queryTrigPoint1 = TrigPoint('AminoAcids', 'M', 10)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 1, 9)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 11)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 20, 9)
        queryTrigPoint2 = TrigPoint('AminoAcids', 'M', 30)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 21, 9)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 30)

        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark1,
            'queryTrigPoint': queryTrigPoint1,
            'subjectLandmark': subjectLandmark1,
            'subjectTrigPoint': subjectTrigPoint1,
        })
        histogram.add(44, {
            'queryLandmark': queryLandmark2,
            'queryTrigPoint': queryTrigPoint2,
            'subjectLandmark': subjectLandmark2,
            'subjectTrigPoint': subjectTrigPoint2,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [AminoAcids])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 20,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'overallScoreAdjustedToBestBinScore': False,
                'queryOffsetsInBins': 31,
                'normalizerQuery': 1.0,
                'normalizerSubject': 1.0,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 30,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOneBinOnePairInBinQueryHasOneHashOutsideMatch(self):
        """
        A match with one bin containing one pair must have an overall score of
        1.0 if the query has an additional pair that is outside the match area
        (because the subject should be used to do the normalisation by length).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        queryTrigPoint = TrigPoint('AminoAcids', 'M', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        subjectTrigPoint = TrigPoint('AminoAcids', 'M', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [AminoAcids])
        query = AARead('id', 300 * 'A' + 'FRRRFRRRFAAAC')
        subject = Subject('id', 30 * 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 10,
                'matchedOffsetCount': 20,
                'matchedQueryOffsetCount': 10,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 10,
                'numeratorQuery': 10,
                'numeratorSubject': 10,
                'overallScoreAdjustedToBestBinScore': False,
                'queryOffsetsInBins': 11,
                'normalizerQuery': 0.5,
                'normalizerSubject': 1.0,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 11,
                'totalOffsetCount': 20,
            },
            analysis)

    def testOneBinOnePairInBinQueryHasTwoPairsOutsideMatch(self):
        """
        A match with one bin containing one pair must have an overall score of
        1.0 if the query has two additional pairs that are outside the match
        area (because the subject should be used to do the normalisation by
        length).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 9)
        queryTrigPoint = TrigPoint('AminoAcids', 'M', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 9)
        subjectTrigPoint = TrigPoint('AminoAcids', 'M', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [AminoAcids])
        query = AARead('id', 'FRRRFRRRF' + ('F' * 200) + 'FRRRFRRRFAAACAAAW')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 21,
                'denominatorSubject': 10,
                'matchedOffsetCount': 20,
                'matchedQueryOffsetCount': 10,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 10,
                'numeratorQuery': 10,
                'numeratorSubject': 10,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 10 / 21,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 11,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 11,
                'totalOffsetCount': 20,
            },
            analysis)

    def testOnePairInBinQuery2Subject1PairOutsideMatch(self):
        """
        A bin containing one pair must have an overall score of 2 / 3 if
        the query has two and the subject one pair that are outside the match
        area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [AminoAcids])
        query = AARead('id', 'FRRRFRRRF' + 'AAACAAAW')
        subject = Subject('id2', 'FRRRFRRRF' + 'AAAC', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 2 / 3},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertAlmostEqual(2 / 3, score)
        self.assertEqual(
            {
                'denominatorQuery': 31,
                'denominatorSubject': 30,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 20 / 31,
                'normalizerSubject': 2 / 3,
                'queryOffsetsInBins': 20,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 20,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOneHashInBinQuery1Subject2PairOutsideMatch(self):
        """
        A bin containing one pair must have an overall score of 2 / 3 if
        the query has one and the subject two pairs that are outside the match
        area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [AminoAcids])
        query = Subject('id2', 'FRRRFRRRF' + 'AAAC', 0)
        subject = AARead('id', 'FRRRFRRRF' + 'AAACAAAW')
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 2 / 3},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertAlmostEqual(2 / 3, score)
        self.assertEqual(
            {
                'denominatorQuery': 30,
                'denominatorSubject': 31,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 2 / 3,
                'normalizerSubject': 20 / 31,
                'queryOffsetsInBins': 20,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 20,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOnePairInBinQueryHasOneUnmatchedPairInsideMatch(self):
        """
        A bin containing one pair where the landmark and trig point do not
        overlap must have the correct overall score if the query has an
        additional non-matching pair that is inside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 50)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 50)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 20 * 'A' + 'FRRRFRRRFAAC')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 42 / 52},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual((21 + 21) / (21 + 21 + 10), score)
        self.assertEqual(
            {
                'denominatorQuery': 31,
                'denominatorSubject': 21,
                'matchedOffsetCount': 42,
                'matchedQueryOffsetCount': 21,
                'matchedRegionScore': 42 / 52,
                'matchedSubjectOffsetCount': 21,
                'numeratorQuery': 31,
                'numeratorSubject': 21,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 1.0,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 51,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 51,
                'totalOffsetCount': 52,
            },
            analysis)

    def testOnePairInBinQueryHasOneUnmatchedPairExactlySpanningMatch(self):
        """
        A bin containing one pair must have a the correct overall score if the
        query has an additional pair that exactly spans the match area but the
        additional pairs' offsets match those of the match (and so do not
        affect the score).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        queryTrigPoint = TrigPoint('Peaks', 'P', 13)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 13)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 'FRRRFRRRFAAAAC')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual((9 + 9) / (9 + 9), score)
        self.assertEqual(
            {
                'denominatorQuery': 10,
                'denominatorSubject': 10,
                'matchedOffsetCount': 20,
                'matchedQueryOffsetCount': 10,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 10,
                'numeratorQuery': 10,
                'numeratorSubject': 10,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 1.0,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 14,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 14,
                'totalOffsetCount': 20,
            },
            analysis)

    def testOnePairInBinQueryHasOneUnmatchedPairExceedingMatch(self):
        """
        A match with one bin containing one pair must have an overall score of
        1.0 if the query has an additional pair that exceeds the match area on
        both sides (because the subject is used for the score normalisation by
        length, the score is not affected).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 2, 5)
        queryTrigPoint = TrigPoint('Peaks', 'P', 5)
        subjectLandmark = Landmark('AlphaHelix', 'A', 2, 5)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 5)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 'FRRRFRRRF' + 20 * 'A' + 'C')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(1.0, score)
        # Note that the landmark in the unmatched hash completely overlaps
        # the alpha helix from offset 2-6 in the query. Because we give
        # priority to AAs that do match, only 5 of the 10 AAs in that
        # unmatched hash get counted as not being matched. For that reason,
        # the denominator of the query is 10, not 15.
        self.assertEqual(
            {
                'denominatorQuery': 10,
                'denominatorSubject': 5,
                'matchedOffsetCount': 10,
                'matchedQueryOffsetCount': 5,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 5,
                'numeratorQuery': 5,
                'numeratorSubject': 5,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 0.5,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 5,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 5,
                'totalOffsetCount': 10,
            },
            analysis)

    def testOnePairInBinQueryHasTwoUnmatchedFeaturesInsideMatch(self):
        """
        A match with one bin containing one pair must have the correct overall
        score if the query has two additional features (making a pair) that are
        inside the match area but which do not overlap the features in the
        pair.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 50)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 50)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [AminoAcids])
        query = AARead('id', 22 * 'A' + 'CAAW')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 42 / 44},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(42 / 44, score)
        self.assertEqual(
            {
                'denominatorQuery': 23,
                'denominatorSubject': 21,
                'matchedOffsetCount': 42,
                'matchedQueryOffsetCount': 21,
                'matchedRegionScore': 42 / 44,
                'matchedSubjectOffsetCount': 21,
                'numeratorQuery': 23,
                'numeratorSubject': 21,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 1.0,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 51,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 51,
                'totalOffsetCount': 44,
            },
            analysis)

    def testOnePairInBinQueryHasOneUnmatchedFeatureOverlappingMatchLeft(self):
        """
        A match with one bin containing one pair must have an overall score of
        1.0 if the query has an additional feature that is only partly inside
        the match area (with the extra feature jutting out on the left).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 'FRRRFRRRFC')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 25,
                'denominatorSubject': 20,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 0.8,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 20,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 20,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOnePairInBinQueryHasOneUnmatchedFeatureOverlappingMatchRight(self):
        """
        A match with one bin containing one pair must have an overall score of
        1.0 if the query has an additional feature that is only partly inside
        the match area (with the extra feature jutting out on the right).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 2, 3)
        queryTrigPoint = TrigPoint('Peaks', 'P', 5)
        subjectLandmark = Landmark('AlphaHelix', 'A', 2, 3)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 5)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 'AAAFRRRFRRRFC')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 11,
                'denominatorSubject': 4,
                'matchedOffsetCount': 8,
                'matchedQueryOffsetCount': 4,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 4,
                'numeratorQuery': 4,
                'numeratorSubject': 4,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 4 / 11,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 4,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 4,
                'totalOffsetCount': 8,
            },
            analysis)

    def testTwoPairs(self):
        """
        A match with one bin containing two pairs must have the correct overall
        score if the query and subject both have an additional feature inside
        their match areas.
        """
        queryLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 10)
        subjectLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 10)

        queryLandmark2 = Landmark('AlphaHelix_pi', 'C', 50, 5)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 60)
        subjectLandmark2 = Landmark('AlphaHelix_pi', 'C', 50, 5)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 60)

        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark1,
            'queryTrigPoint': queryTrigPoint1,
            'subjectLandmark': subjectLandmark1,
            'subjectTrigPoint': subjectTrigPoint1,
        })
        histogram.add(44, {
            'queryLandmark': queryLandmark2,
            'queryTrigPoint': queryTrigPoint2,
            'subjectLandmark': subjectLandmark2,
            'subjectTrigPoint': subjectTrigPoint2,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 20 * 'A' + 'FRRRFRRRFC')
        subject = Subject('id2', 25 * 'A' + 'FRRRFRRRFRRRFAAC', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 20 / 44},
        ]
        sbs = SignificantBinScore(significantBins, query, subject, params)
        score, analysis = sbs.calculateScore()
        matched = (3 + 1) + (5 + 1) + (3 + 1) + (5 + 1)
        total = matched + (9 + 1) + (13 + 1)
        self.assertEqual(matched / total, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 24,
                'matchedOffsetCount': 20,
                'matchedQueryOffsetCount': 10,
                'matchedRegionScore': 20 / 44,
                'matchedSubjectOffsetCount': 10,
                'numeratorQuery': 20,
                'numeratorSubject': 24,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 1.0,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 59,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 59,
                'totalOffsetCount': 44,
            },
            analysis)

    def testThreeSignBinsWithOnePairEach(self):
        """
        A match that consists of three significant bins with one pair in each
        bin must return an overall score of 1.0.
        """
        queryLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 10)
        subjectLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 10)

        queryLandmark2 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 60)
        subjectLandmark2 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 60)

        queryLandmark3 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        queryTrigPoint3 = TrigPoint('Peaks', 'P', 70)
        subjectLandmark3 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        subjectTrigPoint3 = TrigPoint('Peaks', 'P', 70)

        histogram = Histogram(3)
        histogram.add(10, {
            'queryLandmark': queryLandmark1,
            'queryTrigPoint': queryTrigPoint1,
            'subjectLandmark': subjectLandmark1,
            'subjectTrigPoint': subjectTrigPoint1,
        })
        histogram.add(20, {
            'queryLandmark': queryLandmark2,
            'queryTrigPoint': queryTrigPoint2,
            'subjectLandmark': subjectLandmark2,
            'subjectTrigPoint': subjectTrigPoint2,
        })
        histogram.add(30, {
            'queryLandmark': queryLandmark3,
            'queryTrigPoint': queryTrigPoint3,
            'subjectLandmark': subjectLandmark3,
            'subjectTrigPoint': subjectTrigPoint3,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 80 * 'A')
        subject = Subject('id2', 80 * 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
            {'index': 1, 'bin': histogram.bins[1], 'score': 1.0},
            {'index': 2, 'bin': histogram.bins[2], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        overallScore, overallAnalysis = sbs.calculateScore()
        self.maxDiff = None
        self.assertEqual(1.0, overallScore)
        self.assertEqual(
            {
                'denominatorQuery': 11,
                'denominatorSubject': 11,
                'matchedOffsetCount': 22,
                'matchedQueryOffsetCount': 11,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 11,
                'numeratorQuery': 11,
                'numeratorSubject': 11,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 1.0,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 40,
                'score': overallScore,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 40,
                'totalOffsetCount': 22,
            },
            overallAnalysis)

    def testThreeBinsOnePairEachQueryHasOneNonMatchingPairOutsideMatch(self):
        """
        A match that consists of three significant bins with one pair in each
        bin and an additional non-matching pair in the query must return an
        overall score of 1.0.
        """
        queryLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 10)
        subjectLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 10)

        queryLandmark2 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 60)
        subjectLandmark2 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 60)

        queryLandmark3 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        queryTrigPoint3 = TrigPoint('Peaks', 'P', 70)
        subjectLandmark3 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        subjectTrigPoint3 = TrigPoint('Peaks', 'P', 70)

        histogram = Histogram(3)
        histogram.add(10, {
            'queryLandmark': queryLandmark1,
            'queryTrigPoint': queryTrigPoint1,
            'subjectLandmark': subjectLandmark1,
            'subjectTrigPoint': subjectTrigPoint1,
        })
        histogram.add(20, {
            'queryLandmark': queryLandmark2,
            'queryTrigPoint': queryTrigPoint2,
            'subjectLandmark': subjectLandmark2,
            'subjectTrigPoint': subjectTrigPoint2,
        })
        histogram.add(30, {
            'queryLandmark': queryLandmark3,
            'queryTrigPoint': queryTrigPoint3,
            'subjectLandmark': subjectLandmark3,
            'subjectTrigPoint': subjectTrigPoint3,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 80 * 'A' + 'FRRRFRRRFAAAC')
        subject = Subject('id2', 80 * 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
            {'index': 1, 'bin': histogram.bins[1], 'score': 1.0},
            {'index': 2, 'bin': histogram.bins[2], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        overallScore, overallAnalysis = sbs.calculateScore()
        self.maxDiff = None
        self.assertEqual(1.0, overallScore)
        self.assertEqual(
            {
                'denominatorQuery': 21,
                'denominatorSubject': 11,
                'matchedOffsetCount': 22,
                'matchedQueryOffsetCount': 11,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 11,
                'numeratorQuery': 11,
                'numeratorSubject': 11,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 0.5238095238095238,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 40,
                'score': overallScore,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 40,
                'totalOffsetCount': 22,
            },
            overallAnalysis)

    def testThreeBinsOnePairEachQuery2Subject1PairOutsideMatch(self):
        """
        A match that consists of three significant bins with one pair in each
        bin and two additional non-matching pairs in the query and one
        additional non-matching pair in the subject must not return an overall
        score of 1.0.
        """
        queryLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 10)
        subjectLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 10)

        queryLandmark2 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 60)
        subjectLandmark2 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 60)

        queryLandmark3 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        queryTrigPoint3 = TrigPoint('Peaks', 'P', 70)
        subjectLandmark3 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        subjectTrigPoint3 = TrigPoint('Peaks', 'P', 70)

        histogram = Histogram(3)
        histogram.add(10, {
            'queryLandmark': queryLandmark1,
            'queryTrigPoint': queryTrigPoint1,
            'subjectLandmark': subjectLandmark1,
            'subjectTrigPoint': subjectTrigPoint1,
        })
        histogram.add(20, {
            'queryLandmark': queryLandmark2,
            'queryTrigPoint': queryTrigPoint2,
            'subjectLandmark': subjectLandmark2,
            'subjectTrigPoint': subjectTrigPoint2,
        })
        histogram.add(30, {
            'queryLandmark': queryLandmark3,
            'queryTrigPoint': queryTrigPoint3,
            'subjectLandmark': subjectLandmark3,
            'subjectTrigPoint': subjectTrigPoint3,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 80 * 'A' + 'FRRRFRRRF' + 'AAACAAAW')
        subject = Subject('id2', 80 * 'A' + 'FRRRFRRRF' + 'AAAC', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 11 / 21},
            {'index': 1, 'bin': histogram.bins[1], 'score': 11 / 21},
            {'index': 2, 'bin': histogram.bins[2], 'score': 11 / 21},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        overallScore, overallAnalysis = sbs.calculateScore()
        self.maxDiff = None
        self.assertEqual(11 / 21, overallScore)
        self.assertEqual(
            {
                'denominatorQuery': 21,
                'denominatorSubject': 21,
                'matchedOffsetCount': 22,
                'matchedQueryOffsetCount': 11,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 11,
                'numeratorQuery': 11,
                'numeratorSubject': 11,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 11 / 21,
                'normalizerSubject': 11 / 21,
                'queryOffsetsInBins': 40,
                'score': overallScore,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 40,
                'totalOffsetCount': 22,
            },
            overallAnalysis)

    def testThreeSignBinsWithOnePairEachWithOneAdditionalFinder(self):
        """
        A match that consists of three significant bins with one pair in each
        bin and bin with index 2 containing a non-matching pair in the query
        must return the correct overall score.
        """
        queryLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 10)
        subjectLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 10)

        queryLandmark2 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 60)
        subjectLandmark2 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 60)

        queryLandmark3 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        queryTrigPoint3 = TrigPoint('Peaks', 'P', 70)
        subjectLandmark3 = Landmark('AlphaHelix_pi', 'C', 40, 5)
        subjectTrigPoint3 = TrigPoint('Peaks', 'P', 70)

        histogram = Histogram(3)
        histogram.add(10, {
            'queryLandmark': queryLandmark1,
            'queryTrigPoint': queryTrigPoint1,
            'subjectLandmark': subjectLandmark1,
            'subjectTrigPoint': subjectTrigPoint1,
        })
        histogram.add(20, {
            'queryLandmark': queryLandmark2,
            'queryTrigPoint': queryTrigPoint2,
            'subjectLandmark': subjectLandmark2,
            'subjectTrigPoint': subjectTrigPoint2,
        })
        histogram.add(30, {
            'queryLandmark': queryLandmark3,
            'queryTrigPoint': queryTrigPoint3,
            'subjectLandmark': subjectLandmark3,
            'subjectTrigPoint': subjectTrigPoint3,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [Peaks])
        query = AARead('id', 40 * 'A' + 'FRRRFRRRFAAC' + 28 * 'A')
        subject = Subject('id2', 80 * 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 0.0},
            {'index': 1, 'bin': histogram.bins[1], 'score': 0.0},
            {'index': 2, 'bin': histogram.bins[2], 'score': 0.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject, params)
        overallScore, overallAnalysis = sbs.calculateScore()
        self.maxDiff = None
        self.assertEqual(22 / 27, overallScore)
        self.assertEqual(
            {
                'denominatorQuery': 16,
                'denominatorSubject': 11,
                'matchedOffsetCount': 22,
                'matchedQueryOffsetCount': 11,
                'matchedRegionScore': 22 / 27,
                'matchedSubjectOffsetCount': 11,
                'numeratorQuery': 16,
                'numeratorSubject': 11,
                'overallScoreAdjustedToBestBinScore': False,
                'normalizerQuery': 1.0,
                'normalizerSubject': 1.0,
                'queryOffsetsInBins': 40,
                'score': overallScore,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 40,
                'totalOffsetCount': 27,
            },
            overallAnalysis)

    def testCompareEqualSequencesScoreMustBeOne(self):
        """
        If a sequence is compared to itself, the overall score must be 1.0.
        This is a real-life test to check that it actually works.
        Note that this test only passes because we force the overallScore to be
        equal to the bestBinScore if the overall Score is lower than the
        bestBinScore.
        """
        pichninde = AARead('pichninde', 'RLKFGLSYKEQVGGNRELYVGDLNTKLTTRLIEDYS'
                                        'ESLMQNMRYTCLNNEKEFERALLDMKSVVRQSGLAV'
                                        'SMDHSKWGPHMSPVIFAALLKGLEFNLKDGSEVPNA'
                                        'AIVNILLWHIHKMVEVPFNVVEAYMKGFLKRGLGMM'
                                        'DKGGCTIAEEFMFGYFEKGKVPSHISSVLDMGQGIL'
                                        'HNTSDLYGLITEQFINYALELCYGVRFISYTSSDDE'
                                        'IMLSLNEAFKFKDRDELNVDLVLDCMEFHYFLSDKL'
                                        'NKFVSPKTVVGTFASEFKSRFFIWSQEVPLLTKFVA'
                                        'AALH')

        db = DatabaseSpecifier().getDatabaseFromKeywords(
            landmarkNames=[
                'AlphaHelix', 'AlphaHelix_3_10', 'AlphaHelix_pi',
                'AminoAcidsLm', 'BetaStrand', 'BetaTurn', 'Prosite'],
            trigPointNames=['AminoAcids', 'Peaks', 'Troughs'],
            distanceBase=1.01, limitPerLandmark=50, minDistance=1,
            maxDistance=100)
        _, subjectIndex, _ = db.addSubject(pichninde)

        findParams = FindParameters(significanceFraction=0.01,
                                    scoreMethod='FeatureAAScore',
                                    overallScoreMethod='SignificantBinScore')
        result = db.find(pichninde, findParams, storeFullAnalysis=True)
        self.assertEqual(1.0, result.analysis[subjectIndex]['bestBinScore'])
        self.assertEqual(1.0,
                         result.analysis[subjectIndex]['overallScore'])

    def testCompareEqualSequencesScoreMustBeOneWithThreeBins(self):
        """
        If a sequence is compared to itself, the overall score must be 1.0
        when there are three significant bins. This test sequence is a
        subsequence of the Pichninde sequence, chosen to give just a small
        number of significant bins.
        Note that this test only passes because we force the overallScore to be
        equal to the bestBinScore if the overall Score is lower than the
        bestBinScore.
        """
        sequence = AARead('id', 'NTKLTTRLIEDYS')

        db = DatabaseSpecifier().getDatabaseFromKeywords(
            landmarkNames=['Prosite'], trigPointNames=['Troughs'],
            distanceBase=1.0, limitPerLandmark=5, minDistance=1,
            maxDistance=100)
        _, subjectIndex, _ = db.addSubject(sequence)

        findParams = FindParameters(significanceFraction=0.01,
                                    scoreMethod='FeatureAAScore',
                                    overallScoreMethod='SignificantBinScore')
        result = db.find(sequence, findParams, storeFullAnalysis=True)
        self.assertEqual(
            3, len(result.analysis[subjectIndex]['significantBins']))
        self.assertEqual(1.0, result.analysis[subjectIndex]['bestBinScore'])
        self.assertEqual(1.0,
                         result.analysis[subjectIndex]['overallScore'])

    def testScoresMustBeSymmetric(self):
        """
        When comparing two sequences, the scores must be the same, no matter
        which one is used as the query or subject.
        """
        golv = AARead('GOLV', 'RVDIFKKNQHGGLREIYVLDLASRIVQLCLEEISRAVCQELPIEMM'
                              'MHPELKLKKPQEHMYKAAISPESYKSNVSSSNDAKVWNQGHHVAKF'
                              'AQFLCRLLSPEWHGLIVNGLKLWTNKKIALPDGVMNILSRANTPLF'
                              'RNSIHQAVHDSYKGITPMRWLRPGETFMRIESGMMQGILHYTSSLF'
                              'HASLLMMRDSLWRSYSEQLGVKSITTDLVSSDDSSRMTDIFYRDSK'
                              'NFKRGKIFARADHMAIEPLSRCFGIWMSPKSTYCCNGIMEFNSEYF'
                              'FRASLYRPTLKWSYACLG')

        akav = AARead('AKAV', 'VFTYFNKGQKTAKDREIFVGEFEAKMCLYLVERISKERCKLNPDEM'
                              'ISEPGDGKLKKLEDMAEYEIRYTANTLKSMKDKALQEFSKFADDFN'
                              'FKPHSTKIEINADMSKWSAQDVLFKYFWLFALDPALYKPEKERILY'
                              'FLCNYMDKVLVIPDDVMTSILDQRVKREKDIIYEMTNGLKQNWVSI'
                              'KRNWLQGNLNYTSSYLHSCCMNVYKDIIKNVATLLEGDVLVNSMVH'
                              'SDDNHTSITMIQDKFPDDIIIEYCIKLFEKICLSFGNQANMKKTYV'
                              'TNFIKEFVSLFNIYGEPFSVYGRFLLTAVG')

        findParams = FindParameters(significanceFraction=0.01,
                                    scoreMethod='FeatureAAScore',
                                    overallScoreMethod='SignificantBinScore')

        kwds = dict(landmarkNames=['Prosite'], trigPointNames=['Peaks'],
                    distanceBase=1, limitPerLandmark=40, minDistance=1,
                    maxDistance=10000)

        db1 = DatabaseSpecifier().getDatabaseFromKeywords(**kwds)
        _, subjectIndex1, _ = db1.addSubject(golv)
        result1 = db1.find(akav, findParams, storeFullAnalysis=True)

        db2 = DatabaseSpecifier().getDatabaseFromKeywords(**kwds)
        _, subjectIndex2, _ = db2.addSubject(akav)
        result2 = db2.find(golv, findParams, storeFullAnalysis=True)

        self.assertEqual(result1.analysis[subjectIndex1]['bestBinScore'],
                         result2.analysis[subjectIndex2]['bestBinScore'])
        self.assertEqual(result1.analysis[subjectIndex1]['overallScore'],
                         result2.analysis[subjectIndex2]['overallScore'])

    def testPrintAnalysis(self):
        """
        The analysis of a score calculation must print correctly, whether
        we print it using the class name explicitly or the score class that's
        given in the analysis.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [])
        query = AARead('id', 'FRRRFRRRFC')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 1.0},
        ]
        sbs = SignificantBinScore(significantBins, query, subject,
                                  params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(1.0, score)
        self.maxDiff = None
        expected = (
            'Overall score method: SignificantBinScore\n'
            'Overall score: 1.0\n'
            'Total (query+subject) AA offsets in matched pairs in all '
            'bins: 40\n'
            'Subject AA offsets in matched pairs in all bins: 20\n'
            'Query AA offsets in matched pairs in all bins: 20\n'
            'Total (query+subject) AA offsets in hashes in matched '
            'region: 40\n'
            'Matched region score 1.0000 (40 / 40)\n'
            'Query normalizer: 0.8000 (20 / 25)\n'
            'Subject normalizer: 1.0000 (20 / 20)\n'
            'Total query offsets that are in a bin: 20\n'
            'Total subject offsets that are in a bin: 20')
        self.assertEqual(expected, SignificantBinScore.printAnalysis(analysis))
        self.assertEqual(expected,
                         analysis['scoreClass'].printAnalysis(analysis))


class TestGreedySignificantBinScore(TestCase):
    """
    Tests for the light.overall_score.GreedySignificantBinScore class.
    """
    def testOnePairInOneBin(self):
        """
        A match with one bin containing one pair must have an overall score of
        1.0 if the query and subject have no additional (non-matching) hashes.
        This test is already used in the test for the FeatureAAScore and the
        BestBinScore. It's just here to make sure the basics work. Note that
        the bestBinScore that is passed in is too low. This is to make sure the
        GreedySignificantBinScore calculates a score of 1.0.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 110)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 0.0},
        ]
        gsbs = GreedySignificantBinScore(significantBins, query, subject,
                                         params)
        score, analysis = gsbs.calculateScore()
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 20,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'queryOffsetsInBins': 20,
                'normalizerQuery': 1.0,
                'normalizerSubject': 1.0,
                'numberOfBinsConsidered': 1,
                'score': score,
                'scoreClass': GreedySignificantBinScore,
                'subjectOffsetsInBins': 20,
                'totalOffsetCount': 40,
            },
            analysis)

    def testCorrectScoresMustBeCalculated(self):
        """
        When comparing two sequences, the correct score must be calculated.
        """
        eel = AAReadWithX('EelVirusEuropeanX',
                          'SLETLDYINVTVVQDLTDCLTSKGKLPAYLGSKTSETTSILQPWEKETKIP'
                          'VIRRAAKLRAAITWFVEPDSLLAQSILNNIESLTGEDWSASISGFKRTGSA'
                          'LHRFTSARVSAGGFSAQSPARLTRMMATTDTFREIGSDNYDFMFQSLLLFA'
                          'QMTTGEIYKRSPATNFHFHLSCHQCLRKIEEPTLNSDFAYNPIQRSDILDK'
                          'WKPQTTDWSSERKAPEIEEGNWDRLTHQEQSFQVGKSIGFLFGDLTMTKNS'
                          'HAQDSSIFPLSIQYKITAAEFLEGILDGIVKASALSTIHRRNFDHHSKYKS'
                          'TVSGTVDYLIELISESAGFTNLTRNGPLKACLTIPHKIPPSYPLSQSDLGA'
                          )
        ves = AAReadWithX('VesicularStomatitisIndianaVirus',
                          'TSGFNYVSVHCPDGIHDVFSSRGPLPAYLGSKTSESTSILQPWERESKVPL'
                          'IKRATRLRDAISWFVEPDSKLAMTILSNIHSLTGEEWTKRQHGFKRTGSAL'
                          'HRFSTSRMSHGGFASQSTAALTRLMATTDTMRDLGDQNFDFLFQATLLYAQ'
                          'ITTTVARDGWITSCTDHYHIACKSCLRPIEEITLDSSMDYTPPDVSHVLKT'
                          'WRNGEGSWGQEIKQIYPLEGNWKNLAPAEQSYQVGRCIGFLYGDLAYRKST'
                          'HAEDSSLFPLSIQGRIRGRGFLKGLLDGLMRASCCQVIHRRSLAHLKRPAN'
                          'AVYGGLIYLIDKLSVSPPFLSLTRSGPIRDELTIPHKIPTSYPTSNRDMGV'
                          )

        findParams = FindParameters(
            significanceFraction=0.01, scoreMethod='FeatureAAScore',
            overallScoreMethod='GreedySignificantBinScore')

        kwds = dict(landmarkNames=[
            'AlphaHelix', 'AlphaHelix_3_10', 'AlphaHelix_pi', 'AminoAcidsLm',
            'BetaStrand', 'BetaTurn', 'Prosite', 'GOR4AlphaHelix',
            'GOR4BetaStrand'],
            trigPointNames=['Peaks', 'Troughs', 'AminoAcids',
                            'IndividualPeaks', 'IndividualTroughs'],
            distanceBase=1.0, limitPerLandmark=50, minDistance=1,
            maxDistance=100, featureLengthBase=1.01)

        db = DatabaseSpecifier().getDatabaseFromKeywords(**kwds)
        _, subjectIndex, _ = db.addSubject(eel)
        result = db.find(ves, findParams, storeFullAnalysis=True)

        self.assertEqual(0.3023133221365436,
                         result.analysis[subjectIndex]['bestBinScore'])
        self.assertAlmostEqual(0.5842663952304784,
                               result.analysis[subjectIndex]['overallScore'])
        self.maxDiff = None
        analysis = result.analysis[subjectIndex]['overallScoreAnalysis']
        self.assertEqual(
            {
                'denominatorQuery': 301,
                'denominatorSubject': 299,
                'matchedOffsetCount': 196,
                'matchedQueryOffsetCount': 97,
                'matchedRegionScore': 196 / 599,
                'matchedSubjectOffsetCount': 99,
                'numeratorQuery': 301,
                'numeratorSubject': 298,
                'queryOffsetsInBins': 355,
                'normalizerQuery': 1.0,
                'normalizerSubject': 298 / 299,
                'numberOfBinsConsidered': 3,
                'score': 0.5842663952304784,
                'scoreClass': GreedySignificantBinScore,
                'subjectOffsetsInBins': 354,
                'totalOffsetCount': 599,
            },
            analysis)

    def testPrintAnalysis(self):
        """
        The analysis of a score calculation must print correctly, whether
        we print it using the class name explicitly or the score class that's
        given in the analysis.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 10)
        histogram = Histogram(3)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.add(24, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix, AminoAcidsLm], [Peaks])
        query = AARead('id', 'FRRRFRRRFC')
        subject = Subject('id2', 'A', 0)
        significantBins = [
            {'index': 0, 'bin': histogram.bins[0], 'score': 0.0},
        ]
        gsbs = GreedySignificantBinScore(significantBins, query, subject,
                                         params)
        score, analysis = gsbs.calculateScore()
        expected = (
            'Overall score method: GreedySignificantBinScore\n'
            'Overall score: 1.0\n'
            'Total (query+subject) AA offsets in matched pairs in all '
            'bins: 40\n'
            'Subject AA offsets in matched pairs in all bins: 20\n'
            'Query AA offsets in matched pairs in all bins: 20\n'
            'Total (query+subject) AA offsets in hashes in matched '
            'region: 40\n'
            'Matched region score 1.0000 (40 / 40)\n'
            'Query normalizer: 0.8000 (20 / 25)\n'
            'Subject normalizer: 1.0000 (20 / 20)\n'
            'Total query offsets that are in a bin: 20\n'
            'Total subject offsets that are in a bin: 20\n'
            'Number of bins included in the score calculation: 1')
        self.assertEqual(expected,
                         GreedySignificantBinScore.printAnalysis(analysis))
        self.assertEqual(expected,
                         analysis['scoreClass'].printAnalysis(analysis))
