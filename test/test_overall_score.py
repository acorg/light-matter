from unittest import TestCase

from dark.reads import AARead

from light.database import DatabaseSpecifier
from light.features import Landmark, TrigPoint
from light.histogram import Histogram
from light.landmarks import AlphaHelix, AminoAcids as AminoAcidsLm
from light.trig.amino_acids import AminoAcids
from light.trig.peaks import Peaks
from light.overall_score import (BestBinScore, SignificantBinScore,
                                 featureIndicesOutsideOffsets)
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


class TestFeatureIndicesOutsideOffsets(TestCase):
    """
    Tests for the light.overall_score.featureIndicesOutsideOffsets function.
    """
    def testReturnsRightOffsetsForLandmarkInside(self):
        """
        The featureIndicesOutsideOffsets function must return the right offsets
        for a landmark completely inside the offsets specified.
        """
        landmark = Landmark('AlphaHelix', 'A', 10, 10)
        offsets = featureIndicesOutsideOffsets(landmark, set(range(0, 21)))
        self.assertEqual(set(), offsets)

    def testReturnsRightOffsetsForLandmarkPartlyInside(self):
        """
        The featureIndicesOutsideOffsets function must return the right offsets
        for a landmark partly inside the offsets specified.
        """
        landmark = Landmark('AlphaHelix', 'A', 10, 10)
        offsets = featureIndicesOutsideOffsets(landmark, set(range(0, 15)))
        self.assertEqual({15, 16, 17, 18, 19}, offsets)

    def testReturnsRightOffsetsForTrigPointInside(self):
        """
        The featureIndicesOutsideOffsets function must return the right offsets
        for a trigPoint completely inside the offsets specified.
        """
        trigPoint = TrigPoint('Peaks', 'P', 10)
        offsets = featureIndicesOutsideOffsets(trigPoint, set(range(0, 5)))
        self.assertEqual({10}, offsets)


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
        significantBinScore = SignificantBinScore(histogram, [], query,
                                                  subject, params)
        score, analysis = significantBinScore.calculateScore()
        self.assertIs(None, score)
        self.assertEqual({}, analysis)

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
        sbs = SignificantBinScore(histogram, [], query, subject, params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(None, score)
        self.assertEqual({}, analysis)

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
        sbs = SignificantBinScore(histogram, [], query, subject, params)
        score, analysis = sbs.calculateScore()
        self.assertEqual(None, score)
        self.assertEqual({}, analysis)

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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'queryOffsetsInBins': 20,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'queryOffsetsInBins': 31,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'queryOffsetsInBins': 11,
                'normaliserQuery': 0.5,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 10 / 21,
                'normaliserSubject': 1.0,
                'queryOffsetsInBins': 11,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 11,
                'totalOffsetCount': 20,
            },
            analysis)

    def testOnePairInBinQuery2Subject1PairOutsideMatch(self):
        """
        A bin containing one pair must not have an overall score of 2 / 3 if
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 20 / 31,
                'normaliserSubject': 2 / 3,
                'queryOffsetsInBins': 20,
                'score': score,
                'scoreClass': SignificantBinScore,
                'subjectOffsetsInBins': 20,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOneHashInBinQuery1Subject2PairOutsideMatch(self):
        """
        A bin containing one pair must not have an overall score of 2 / 3 if
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 2 / 3,
                'normaliserSubject': 20 / 31,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 0.5,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 0.8,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 4 / 11,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
                                  params)
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
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}, {'index': 1}, {'index': 2}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
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
        significantBins = [{'index': 0}, {'index': 1}, {'index': 2}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 0.5238095238095238,
                'normaliserSubject': 1.0,
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
        additional non-matching pair in the query must not return an overall
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
        significantBins = [{'index': 0}, {'index': 1}, {'index': 2}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
                'normaliserQuery': 11 / 21,
                'normaliserSubject': 11 / 21,
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
        significantBins = [{'index': 0}, {'index': 1}, {'index': 2}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
                                  params)
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
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
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
        self.assertEqual(0.5060240963855421,
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
        significantBins = [{'index': 0}]
        sbs = SignificantBinScore(histogram, significantBins, query, subject,
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
