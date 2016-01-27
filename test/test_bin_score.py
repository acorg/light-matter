from __future__ import division

import six
import warnings
from unittest import TestCase

from dark.reads import AARead

from light.database import DatabaseSpecifier
from light.features import Landmark, TrigPoint
from light.parameters import Parameters, FindParameters
from light.subject import Subject
from light.bin_score import (
    NoneScore, MinHashesScore, FeatureMatchingScore, FeatureAAScore,
    histogramBinFeatures, featureInRange, getHashFeatures,
    weightedHistogramBinFeatures, getWeightedOffsets, WeightedFeatureAAScore)
from light.histogram import Histogram
from light.landmarks import AlphaHelix, AminoAcids as AminoAcidsLm
from light.trig import Peaks, AminoAcids

DEFAULT_WEIGHTS = {
    'AlphaHelix': 1,
    'AlphaHelix_3_10': 1,
    'AlphaHelix_pi': 1,
    'BetaStrand': 1,
    'BetaTurn': 1,
    'AminoAcidsLm': 1,
    'GOR4AlphaHelix': 1,
    'GOR4BetaStrand': 1,
    'GOR4Coil': 1,
    'Prosite': 1,
    'Peaks': 1,
    'Troughs': 1,
    'AminoAcids': 1,
    'IndividualPeaks': 1,
    'IndividualTroughs': 1,
}

TEST_WEIGHTS = {
    'AlphaHelix': 1.5,
    'AlphaHelix_3_10': 1,
    'AlphaHelix_pi': 1,
    'BetaStrand': 1,
    'BetaTurn': 1,
    'AminoAcidsLm': 1,
    'GOR4AlphaHelix': 1,
    'GOR4BetaStrand': 1,
    'GOR4Coil': 1,
    'Prosite': 1,
    'Peaks': 1,
    'Troughs': 1,
    'AminoAcids': 1,
    'IndividualPeaks': 1,
    'IndividualTroughs': 1,
}


class TestNoneScore(TestCase):
    """
    Tests for the light.bin_score.NoneScore class.
    """
    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of C{None}.
        """
        histogram = Histogram()
        histogram.finalize()
        ns = NoneScore()
        score, analysis = ns.calculateScore(0)
        self.assertIs(None, score)
        self.assertEqual(
            {
                'score': None,
                'scoreClass': NoneScore,
            },
            analysis)

    def testOneHashInBin(self):
        """
        A bin containing one hash must have a score of C{None}.
        """
        # Note that this test is almost identical to a FeatureMatchingScore
        # test below (with the same test name), but will result in a score
        # of None.
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
        ns = NoneScore()
        score, analysis = ns.calculateScore(0)
        self.assertIs(None, score)
        self.assertEqual(
            {
                'score': None,
                'scoreClass': NoneScore,
            },
            analysis)


class TestMinHashesScore(TestCase):
    """
    Tests for the light.bin_score.MinHashesScore class.
    """
    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of 0.0
        """
        histogram = Histogram()
        histogram.finalize()
        mhs = MinHashesScore(histogram, 0)
        score, analysis = mhs.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'binCount': 0,
                'minHashCount': 0,
                'score': score,
                'scoreClass': MinHashesScore,
            },
            analysis)

    def testZeroMinHashCount(self):
        """
        If the minimum hash count for the query/subject is zero, the score
        for any bin must be 0.0
        """
        # Use a histogram with one bin, and put something in it, to make sure
        # we're not testing an empty histogram bin.
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 0)
        score, analysis = mhs.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'binCount': 1,
                'minHashCount': 0,
                'score': score,
                'scoreClass': MinHashesScore,
            },
            analysis)

    def testBinWithOneElementMinHashCountOne(self):
        """
        If the minimum hash count for the query/subject is one and the bin has
        one thing in it, the score for the bin must be 1.0
        """
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 1)
        score, analysis = mhs.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'binCount': 1,
                'minHashCount': 1,
                'score': score,
                'scoreClass': MinHashesScore,
            },
            analysis)

    def testBinWithOneElementMinHashCountTwo(self):
        """
        If the minimum hash count for the query/subject is two and the bin has
        one thing in it, the score for the bin must be 0.5
        """
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 2)
        score, analysis = mhs.calculateScore(0)
        self.assertEqual(0.5, score)
        self.assertEqual(
            {
                'binCount': 1,
                'minHashCount': 2,
                'score': score,
                'scoreClass': MinHashesScore,
            },
            analysis)

    def testBinWithTwoElementsMinHashCountOne(self):
        """
        If the minimum hash count for the query/subject is one and the bin has
        two things in it, a warning will be issued to say that the bin has more
        things in it than should be possible (given the minimum hash count of
        the query/subject).
        """
        histogram = Histogram(1)
        histogram.add(444)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 1)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            score, analysis = mhs.calculateScore(0)
            self.assertEqual(1, len(w))
            self.assertTrue(issubclass(w[0].category, RuntimeWarning))
            error = ('Bin contains 2 deltas for a query/subject pair with a '
                     'minimum hash count of only 1.')
            self.assertIn(error, str(w[0].message))
        self.assertEqual(
            {
                'binCount': 2,
                'minHashCount': 1,
                'score': score,
                'scoreClass': MinHashesScore,
            },
            analysis)

    def testPrintAnalysis(self):
        """
        The analysis of a score calculation must print correctly, whether
        we print it using the class name explicitly or the score class that's
        given in the analysis.
        """
        histogram = Histogram(1)
        histogram.add(444)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 1)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            score, analysis = mhs.calculateScore(0)
            self.assertEqual(1, len(w))
            self.assertTrue(issubclass(w[0].category, RuntimeWarning))
            error = ('Bin contains 2 deltas for a query/subject pair with a '
                     'minimum hash count of only 1.')
            self.assertIn(error, str(w[0].message))

        expected = (
            'Score method: MinHashesScore\n'
            'Minimum hash count: 1\n'
            'Bin count: 2\n'
            'Score: %.4f' % score)

        self.assertEqual(expected, MinHashesScore.printAnalysis(analysis))
        self.assertEqual(expected,
                         analysis['scoreClass'].printAnalysis(analysis))


class TestHistogramBinFeatures(TestCase):
    """
    Tests for the light.bin_score.histogramBinFeatures function.
    """
    def testInvalidQueryOrSubjectSpecifier(self):
        """
        If an invalid value for queryOrSubject ('xxx') is passed to
        histogramBinFeatures, it must raise a KeyError (when looking
        in the bin item dictionary for 'xxxLandmark').
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 111)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        six.assertRaisesRegex(self, KeyError, 'xxxLandmark',
                              histogramBinFeatures, histogram[0], 'xxx')

    def testEmptyBin(self):
        """
        If a histogram bin is empty, then when asked to extract its features
        and offsets, histogramBinFeatures must return two empty sets.
        """
        histogram = Histogram()
        histogram.finalize()
        features, offsets = histogramBinFeatures(histogram[0], 'query')
        self.assertEqual(set(), features)
        self.assertEqual(set(), offsets)

    def testOneFeatureQuery(self):
        """
        If a histogram bin has just one feature, histogramBinFeatures must
        return the details of that feature in the query.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 111)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        features, offsets = histogramBinFeatures(histogram[0], 'query')
        self.assertEqual(set([queryLandmark, queryTrigPoint]), features)
        self.assertEqual(set(range(100, 120)), offsets)

    def testOneFeatureSubject(self):
        """
        If a histogram bin has just one feature, histogramBinFeatures must
        return the details of that feature in the subject.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 111)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        features, offsets = histogramBinFeatures(histogram[0], 'subject')
        self.assertEqual(set([subjectLandmark, subjectTrigPoint]), features)
        self.assertEqual(set(range(101, 121)), offsets)

    def testOneFeatureTwoLocations(self):
        """
        If a histogram bin has one feature that appears in two places,
        histogramBinFeatures must return the details of that feature.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 111)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 210)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 201, 20)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 211)

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
        features, offsets = histogramBinFeatures(histogram[0], 'query')
        self.assertEqual(set([queryLandmark1, queryTrigPoint1,
                              queryLandmark2, queryTrigPoint2]),
                         features)
        expectedOffsets = set(range(100, 120))
        expectedOffsets.update(range(200, 220))
        self.assertEqual(expectedOffsets, offsets)


class TestFeatureInRange(TestCase):
    """
    Tests for the light.bin_score.featureInRange function.
    """
    def testNoMatchWithNoneMinOffset(self):
        """
        If the min offset is passed as None, featureInRange must return False.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        self.assertFalse(featureInRange(landmark, None, 2000))

    def testNoMatchLeft(self):
        """
        If the feature is to the left of the min offset, featureInRange
        must return False.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        self.assertFalse(featureInRange(landmark, 1000, 2000))

    def testNoMatchRight(self):
        """
        If the feature is to the right of the min offset, featureInRange
        must return False.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        self.assertFalse(featureInRange(landmark, 0, 10))

    def testNoMatchLeftOverlap(self):
        """
        If the feature extends partly to the left of the min offset,
        featureInRange must return False.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        self.assertFalse(featureInRange(landmark, 110, 200))

    def testNoMatchRightOverlap(self):
        """
        If the feature extends partly to the right of the min offset,
        featureInRange must return False.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        self.assertFalse(featureInRange(landmark, 0, 110))

    def testExactOverlap(self):
        """
        If the feature exactly matches the allowed min/max offsets,
        featureInRange must return True.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        self.assertTrue(featureInRange(landmark, 100, 119))

    def testLeftAligned(self):
        """
        If the feature starts at the min offset and ends within the allowed
        maximum offset, featureInRange must return True.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        self.assertTrue(featureInRange(landmark, 100, 200))

    def testRightAligned(self):
        """
        If the feature ends at the max offset and starts after the allowed
        minimum offset, featureInRange must return True.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        self.assertTrue(featureInRange(landmark, 10, 119))

    def testFullyContained(self):
        """
        If the feature is completely within the min/max allowed offset,
        featureInRange must return True.
        """
        landmark = Landmark('AlphaHelix', 'A', 100, 20)
        self.assertTrue(featureInRange(landmark, 10, 200))


class TestGetHashFeatures(TestCase):
    """
    Tests for the light.bin_score.getHashFeatures function.
    """
    def testNoHashes(self):
        """
        If no hashes are passed in, the returned set must be empty.
        """
        hashes = {}
        result = getHashFeatures(hashes)
        self.assertEqual(set(), result)

    def testOneHashOneLocation(self):
        """
        If one hash in one location is passed in, two features must be
        returned.
        """
        helixAt0 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 0, 9, 2)
        peakAt10 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 10)

        hashes = {
            'A2:P:15': [[helixAt0, peakAt10]],
        }

        result = getHashFeatures(hashes)
        self.assertEqual(set([helixAt0, peakAt10]), result)

    def testOneHashTwoLocations(self):
        """
        If a hash which occurs in two locations is passed in, four features
        with the right offsets must be returned.
        """
        helixAt0 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 0, 9, 2)
        peakAt10 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 10)
        helixAt30 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 30, 9, 2)
        peakAt40 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 40)

        hashes = {
            'A2:A2:15': [[helixAt0, peakAt10], [helixAt30, peakAt40]]
        }

        result = getHashFeatures(hashes)
        self.assertEqual(set([helixAt0, peakAt10, helixAt30, peakAt40]),
                         result)

    def testTwoHashesTwoLocations(self):
        """
        If two hashes in two different locations are passed in, the four
        correct features must be returned.
        """
        helixAt0 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 0, 9, 2)
        peakAt10 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 10)
        helixAt15 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 15, 9, 2)
        peakAt13 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 13)

        hashes = {
            'A2:A2:15': [[helixAt0, peakAt10]],
            'A2:P:-2': [[helixAt15, peakAt13]],
        }

        result = getHashFeatures(hashes)
        self.assertEqual(set([helixAt0, peakAt10, helixAt15, peakAt13]),
                         result)

    def testTwoHashesOneOfThemInTwoLocations(self):
        """
        If two hashes, one of which occurs in two locations are passed in, six
        features with the correct offsets must be returned.
        """
        helixAt0 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 0, 9, 2)
        peakAt10 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 10)
        helixAt30 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 30, 9, 2)
        peakAt40 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 40)
        helixAt15 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 15, 9, 2)
        peakAt13 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 13)

        hashes = {
            'A2:A2:15': [[helixAt0, peakAt10], [helixAt30, peakAt40]],
            'A2:P:-2': [[helixAt15, peakAt13]],
        }

        result = getHashFeatures(hashes)
        self.assertEqual(set([helixAt0, peakAt10, helixAt30, peakAt40,
                         helixAt15, peakAt13]),
                         result)


class TestFeatureMatchingScore(TestCase):
    """
    Tests for the light.bin_score.FeatureMatchingScore class.
    """
    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of 0.0 if the query and
        subject both have no features.
        """
        histogram = Histogram()
        histogram.finalize()
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'matchScore': 0.0,
                'maxQueryOffset': None,
                'maxSubjectOffset': None,
                'minQueryOffset': None,
                'minSubjectOffset': None,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testEmptyBinQueryHasOneFeature(self):
        """
        A bin containing no hashes must have a score of zero, even if the query
        has a feature. There is no match region, so no unmatched features can
        fall inside it.
        """
        histogram = Histogram(1)
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'matchScore': 0.0,
                'maxQueryOffset': None,
                'maxSubjectOffset': None,
                'minQueryOffset': None,
                'minSubjectOffset': None,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBin(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four) if the query and subject have no
        additional (non-matching) features. The reward is multiplied by four
        because both the landmark and the trig point (2) appear in both the
        query and the subject (2) and 2 x 2 = 4.
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
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinOccurringInTwoPlaces(self):
        """
        A bin containing one hash that occurs in two place must have a score
        that is the feature match reward multiplied by eight) if the query and
        subject have no additional (non-matching) features. The reward is
        multiplied by eight because both the landmark and the trig point (2)
        appear twice (2) in both the query and the subject (2) and
        2 x 2 x 2 = 8.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 111)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 210)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 201, 20)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 211)

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
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(8 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         score)
        self.assertEqual(
            {
                'matchScore': 8 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 219,
                'maxSubjectOffset': 220,
                'minQueryOffset': 100,
                'minSubjectOffset': 101,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryHasOneFeatureOutsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four if the query has an additional feature
        that is outside the match area.
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
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryHasTwoFeaturesOutsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four if the query has two additional
        features that are outside the match area.
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
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF' + ('F' * 200) + 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryAndSubjectHaveOneFeaturesOutsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four if the query and the subject each have
        one additional feature that are both outside the match area.
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
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', ('F' * 200) + 'FRRRFRRRF', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureInsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus the feature mismatch score if the
        query has an additional feature that is inside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE +
            FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 19,
                'maxSubjectOffset': 19,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'mismatchScore': FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureExactlySpanningMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus the feature mismatch score if the
        query has an additional feature that exactly spans the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        queryTrigPoint = TrigPoint('Peaks', 'P', 5)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 5)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE +
            FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 8,
                'maxSubjectOffset': 8,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'mismatchScore': FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureExceedingMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four if the query has an additional feature that
        exceeds the match area on both sides.
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
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
            score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 6,
                'maxSubjectOffset': 6,
                'minQueryOffset': 2,
                'minSubjectOffset': 2,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryHasTwoUnmatchedFeatureInsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus twice the feature mismatch score if
        the query has two additional features that are inside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRFAFRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE +
            2 * FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 19,
                'maxSubjectOffset': 19,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'mismatchScore': (
                    2 * FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE),
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchLeft(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four if the query has an additional feature that
        is only partly inside the match area (with the extra feature jutting
        out on the left).
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
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
            score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 24,
                'maxSubjectOffset': 24,
                'minQueryOffset': 5,
                'minSubjectOffset': 5,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchRight(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four if the query has an additional feature that
        is only partly inside the match area (with the extra feature jutting
        out on the right).
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
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'AAAFRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
            score)
        self.assertEqual(
            {
                'matchScore': 4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                'maxQueryOffset': 5,
                'maxSubjectOffset': 5,
                'minQueryOffset': 2,
                'minSubjectOffset': 2,
                'mismatchScore': 0.0,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureInMatchNonDefault(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus the feature mismatch score if the
        query has an additional feature that is inside the match area,
        including when non-default values are used for the feature match and
        mismatch scores.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        findParams = FindParameters(featureMatchScore=3.1,
                                    featureMismatchScore=-1.2)
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params,
                                   findParams)
        score, analysis = fms.calculateScore(0)
        self.assertEqual(4 * 3.1 - 1.2, score)
        self.assertEqual(
            {
                'matchScore': 4 * 3.1,
                'maxQueryOffset': 19,
                'maxSubjectOffset': 19,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'mismatchScore': -1.2,
                'score': score,
                'scoreClass': FeatureMatchingScore,
            },
            analysis)

    def testPrintAnalysis(self):
        """
        The analysis of a score calculation must print correctly, whether
        we print it using the class name explicitly or the score class that's
        given in the analysis.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 10)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        findParams = FindParameters(featureMatchScore=3.1,
                                    featureMismatchScore=-1.2)
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        fms = FeatureMatchingScore(histogram, query, subject, params,
                                   findParams)
        score, analysis = fms.calculateScore(0)
        expected = (
            'Score method: FeatureMatchingScore\n'
            'Matched offset range in query: 0 to 19\n'
            'Matched offset range in subject: 0 to 19\n'
            'Match score: 12.4000\n'
            'Mismatch score: -1.2000\n'
            'Score: 11.2000')

        self.assertEqual(expected,
                         FeatureMatchingScore.printAnalysis(analysis))
        self.assertEqual(expected,
                         analysis['scoreClass'].printAnalysis(analysis))


class TestFeatureAAScore(TestCase):
    """
    Tests for the light.bin_score.FeatureAAScore class.
    """

    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of 0.0 if the query and
        subject both have no features.
        """
        histogram = Histogram()
        histogram.finalize()
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 0,
                'denominatorSubject': 0,
                'matchedOffsetCount': 0,
                'matchedQueryOffsetCount': 0,
                'matchedRegionScore': 0.0,
                'matchedSubjectOffsetCount': 0,
                'maxQueryOffset': None,
                'maxSubjectOffset': None,
                'minQueryOffset': None,
                'minSubjectOffset': None,
                'numeratorQuery': 0,
                'numeratorSubject': 0,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 0,
            },
            analysis)

    def testEmptyBinQueryHasOneFeature(self):
        """
        A bin containing no hashes must have a score of zero, even if the query
        has a feature (but no hashes). There is no match region, so that part
        of the score is zero which causes the overall score to be zero.
        """
        histogram = Histogram(1)
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 0,
                'denominatorSubject': 0,
                'matchedOffsetCount': 0,
                'matchedQueryOffsetCount': 0,
                'matchedRegionScore': 0.0,
                'matchedSubjectOffsetCount': 0,
                'maxQueryOffset': None,
                'maxSubjectOffset': None,
                'minQueryOffset': None,
                'minSubjectOffset': None,
                'numeratorQuery': 0,
                'numeratorSubject': 0,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 0,
            },
            analysis)

    def testEmptyBinQueryAndSubjectHaveOneFeature(self):
        """
        A bin containing no hashes must have a score of zero, when the query
        and subject both have one feature (but no hashes). There is no match
        region, so that part of the score is zero which causes the overall
        score to be zero.
        """
        histogram = Histogram(1)
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'AAAAAAAAAAAAAAFRRRFRRRF', 0)
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 0,
                'denominatorSubject': 0,
                'matchedOffsetCount': 0,
                'matchedQueryOffsetCount': 0,
                'matchedRegionScore': 0.0,
                'matchedSubjectOffsetCount': 0,
                'maxQueryOffset': None,
                'maxSubjectOffset': None,
                'minQueryOffset': None,
                'minSubjectOffset': None,
                'numeratorQuery': 0,
                'numeratorSubject': 0,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 0,
            },
            analysis)

    def testOneHashInBin(self):
        """
        A bin containing one hash must have a score of 1.0 if the query and
        subject have no additional (non-matching) hashes.
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 20,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOneHashInBinOccurringInTwoPlaces(self):
        """
        A bin containing one hash that occurs in two places must have a score
        of 1.0 if the query and subject have no additional (non-matching)
        hashes.
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 20,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 30,
                'maxSubjectOffset': 30,
                'minQueryOffset': 0,
                'minSubjectOffset': 1,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOneHashInBinQueryHasOneHashOutsideMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has an
        additional hash that is outside the match area (because the subject
        should be used to do the normalisation by length).
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 10,
                'matchedOffsetCount': 20,
                'matchedQueryOffsetCount': 10,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 10,
                'maxQueryOffset': 10,
                'maxSubjectOffset': 10,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'numeratorQuery': 10,
                'numeratorSubject': 10,
                'normaliserQuery': 0.5,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 20,
            },
            analysis)

    def testOneHashInBinQueryHasTwoHashesOutsideMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has two
        additional hashes that are outside the match area (because the subject
        should be used to do the normalisation by length).
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 21,
                'denominatorSubject': 10,
                'matchedOffsetCount': 20,
                'matchedQueryOffsetCount': 10,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 10,
                'maxQueryOffset': 110,
                'maxSubjectOffset': 110,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'numeratorQuery': 10,
                'numeratorSubject': 10,
                'normaliserQuery': 10 / 21,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 20,
            },
            analysis)

    def testOneHashInBinQuery2Subject1HashOutsideMatch(self):
        """
        A bin containing one hash must not have a score of 1.0 if the query has
        two and the subject one hash that are outside the match area.
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertAlmostEqual(2 / 3, score)
        self.assertEqual(
            {
                'denominatorQuery': 31,
                'denominatorSubject': 30,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 20 / 31,
                'normaliserSubject': 2 / 3,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOneHashInBinQuery1Subject2HashOutsideMatch(self):
        """
        A bin containing one hash must not have a score of 1.0 if the query has
        one and the subject two hashes that are outside the match area.
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertAlmostEqual(2 / 3, score)
        self.assertEqual(
            {
                'denominatorQuery': 30,
                'denominatorSubject': 31,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 2 / 3,
                'normaliserSubject': 20 / 31,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedHashInsideMatch(self):
        """
        A bin containing one hash where the landmark and trig point do not
        overlap must have the correct score if the query has an additional
        non-matching hash that is inside the match area.
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual((21 + 21) / (21 + 21 + 10), score)
        self.assertEqual(
            {
                'denominatorQuery': 31,
                'denominatorSubject': 21,
                'matchedOffsetCount': 42,
                'matchedQueryOffsetCount': 21,
                'matchedRegionScore': 42 / 52,
                'matchedSubjectOffsetCount': 21,
                'maxQueryOffset': 50,
                'maxSubjectOffset': 50,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'numeratorQuery': 31,
                'numeratorSubject': 21,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 52,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedHashExactlySpanningMatch(self):
        """
        A bin containing one hash must have a the correct score if the query
        has an additional hash that exactly spans the match area but the
        additional hashes' offsets match those of the match (and so do not
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual((9 + 9) / (9 + 9), score)
        self.assertEqual(
            {
                'denominatorQuery': 10,
                'denominatorSubject': 10,
                'matchedOffsetCount': 20,
                'matchedQueryOffsetCount': 10,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 10,
                'maxQueryOffset': 13,
                'maxSubjectOffset': 13,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'numeratorQuery': 10,
                'numeratorSubject': 10,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 20,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedHashExceedingMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has
        an additional hash that exceeds the match area on both sides (because
        the subject is used for the score normalisation by length, the score
        is not affected).
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
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
                'maxQueryOffset': 6,
                'maxSubjectOffset': 6,
                'minQueryOffset': 2,
                'minSubjectOffset': 2,
                'numeratorQuery': 5,
                'numeratorSubject': 5,
                'normaliserQuery': 0.5,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 10,
            },
            analysis)

    def testOneHashInBinQueryHasTwoUnmatchedFeaturesInsideMatch(self):
        """
        A bin containing one hash must have the correct score if the query has
        two additional features (making a hash) that are inside the match area
        but which do not overlap the features in the hash.
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(42 / 44, score)
        self.assertEqual(
            {
                'denominatorQuery': 23,
                'denominatorSubject': 21,
                'matchedOffsetCount': 42,
                'matchedQueryOffsetCount': 21,
                'matchedRegionScore': 42 / 44,
                'matchedSubjectOffsetCount': 21,
                'maxQueryOffset': 50,
                'maxSubjectOffset': 50,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'numeratorQuery': 23,
                'numeratorSubject': 21,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 44,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchLeft(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has an
        additional feature that is only partly inside the match area (with the
        extra feature jutting out on the left).
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 25,
                'denominatorSubject': 20,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 24,
                'maxSubjectOffset': 24,
                'minQueryOffset': 5,
                'minSubjectOffset': 5,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 0.8,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchRight(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has an
        additional feature that is only partly inside the match area (with the
        extra feature jutting out on the right).
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 11,
                'denominatorSubject': 4,
                'matchedOffsetCount': 8,
                'matchedQueryOffsetCount': 4,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 4,
                'maxQueryOffset': 5,
                'maxSubjectOffset': 5,
                'minQueryOffset': 2,
                'minSubjectOffset': 2,
                'numeratorQuery': 4,
                'numeratorSubject': 4,
                'normaliserQuery': 4 / 11,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 8,
            },
            analysis)

    def testTwoHashes(self):
        """
        A bin containing two hashes must have the correct score if the query
        and subject both have an additional feature inside their match areas.
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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
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
                'maxQueryOffset': 60,
                'maxSubjectOffset': 60,
                'minQueryOffset': 2,
                'minSubjectOffset': 2,
                'numeratorQuery': 20,
                'numeratorSubject': 24,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 44,
            },
            analysis)

    def testCompareEqualSequencesScoreMustBeOne(self):
        """
        If a sequence is compared to itself, the score must be 1.0. See
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
                                    scoreMethod='FeatureAAScore')
        result = db.find(pichninde, findParams, storeFullAnalysis=True)
        self.assertEqual(1.0, result.analysis[subjectIndex]['bestBinScore'])

        scoreAnalysis = result.analysis[
            subjectIndex]['significantBins'][0]['scoreAnalysis']
        self.maxDiff = None
        self.assertEqual(
            {
                'denominatorQuery': 210,
                'denominatorSubject': 210,
                'matchedOffsetCount': 420,
                'matchedQueryOffsetCount': 210,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 210,
                'maxQueryOffset': 290,
                'maxSubjectOffset': 290,
                'minQueryOffset': 1,
                'minSubjectOffset': 1,
                'numeratorQuery': 210,
                'numeratorSubject': 210,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': 1.0,
                'scoreClass': FeatureAAScore,
                'totalOffsetCount': 420,
            },
            scoreAnalysis)

    def testScoresMustBeSymmetric(self):
        """
        When comparing two sequences, the scores must be the same, no matter
        which one is used as the query or subject.

        This was a (formerly) failing test built during the resolution of
        https://github.com/acorg/light-matter/issues/341 based on two of the
        sequences received from Sandra Junglen on March 13, 2015.
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

        db1 = DatabaseSpecifier().getDatabaseFromKeywords(**kwds)
        _, subjectIndex1, _ = db1.addSubject(golv)
        result1 = db1.find(akav, findParams, storeFullAnalysis=True)

        db2 = DatabaseSpecifier().getDatabaseFromKeywords(**kwds)
        _, subjectIndex2, _ = db2.addSubject(akav)
        result2 = db2.find(golv, findParams, storeFullAnalysis=True)

        self.assertEqual(result1.analysis[subjectIndex1]['bestBinScore'],
                         result2.analysis[subjectIndex2]['bestBinScore'])

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
        faas = FeatureAAScore(histogram, query, subject, params)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)

        expected = (
            'Score method: FeatureAAScore\n'
            'Matched offset range in query: 5 to 24\n'
            'Matched offset range in subject: 5 to 24\n'
            'Total (query+subject) AA offsets in matched hashes: 40\n'
            'Subject AA offsets in matched hashes: 20\n'
            'Query AA offsets in matched hashes: 20\n'
            'Total (query+subject) AA offsets in hashes in matched '
            'region: 40\n'
            'Matched region score 1.0000 (40 / 40)\n'
            'Query normalizer: 0.8000 (20 / 25)\n'
            'Subject normalizer: 1.0000 (20 / 20)\n'
            'Score: 1.0000')
        self.assertEqual(expected, FeatureAAScore.printAnalysis(analysis))
        self.assertEqual(expected,
                         analysis['scoreClass'].printAnalysis(analysis))


class TestWeightedHistogramBinFeatures(TestCase):
    """
    Tests for the light.bin_score.weightedHistogramBinFeatures function.
    """
    def testWeightedInvalidQueryOrSubjectSpecifier(self):
        """
        If an invalid value for queryOrSubject ('xxx') is passed to
        weightedHistogramBinFeatures, it must raise a KeyError (when looking
        in the bin item dictionary for 'xxxLandmark').
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 111)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        six.assertRaisesRegex(self, KeyError, 'xxxLandmark',
                              weightedHistogramBinFeatures, histogram[0],
                              'xxx', DEFAULT_WEIGHTS)

    def testWeightedEmptyBin(self):
        """
        If a histogram bin is empty, then when asked to extract its features
        and offsets, weightedHistogramBinFeatures must return an empty set and
        an empty dict.
        """
        histogram = Histogram()
        histogram.finalize()
        features, offsets = weightedHistogramBinFeatures(histogram[0], 'query',
                                                         DEFAULT_WEIGHTS)
        self.assertEqual(set(), features)
        self.assertEqual({}, offsets)

    def testWeightedOneFeatureQuery(self):
        """
        If a histogram bin has just one feature, weightedHistogramBinFeatures
        must return the details of that feature in the query.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 5)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 101, 5)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 111)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        features, offsets = weightedHistogramBinFeatures(histogram[0], 'query',
                                                         DEFAULT_WEIGHTS)
        self.assertEqual(set([queryLandmark, queryTrigPoint]), features)
        self.assertEqual({100: [1], 101: [1], 102: [1],
                          103: [1], 104: [1], 110: [1]}, offsets)

    def testWeightedOneFeatureSubject(self):
        """
        If a histogram bin has just one feature, weightedHistogramBinFeatures
        must return the details of that feature in the subject.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 5)
        queryTrigPoint = TrigPoint('Peaks', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 101, 5)
        subjectTrigPoint = TrigPoint('Peaks', 'P', 111)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        features, offsets = weightedHistogramBinFeatures(histogram[0],
                                                         'subject',
                                                         DEFAULT_WEIGHTS)
        self.assertEqual(set([subjectLandmark, subjectTrigPoint]), features)
        self.assertEqual({101: [1], 102: [1], 103: [1],
                          104: [1], 105: [1], 111: [1]}, offsets)

    def testWeightedOneFeatureTwoLocations(self):
        """
        If a histogram bin has one feature that appears in two places,
        weightedHistogramBinFeatures must return the details of that feature.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 100, 5)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 101, 5)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 111)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 200, 5)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 210)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 201, 5)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 211)

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
        features, offsets = weightedHistogramBinFeatures(histogram[0], 'query',
                                                         DEFAULT_WEIGHTS)
        self.assertEqual(set([queryLandmark1, queryTrigPoint1,
                              queryLandmark2, queryTrigPoint2]),
                         features)
        self.assertEqual({100: [1], 101: [1], 102: [1],
                          103: [1], 104: [1], 110: [1],
                          200: [1], 201: [1], 202: [1],
                          203: [1], 204: [1], 210: [1]}, offsets)

    def testWeightedOneFeatureTwoLocationsNonDefaultWeights(self):
        """
        If non-default weights are used, the right features, offsets and their
        associated weights must be returned.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 100, 5)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 101, 5)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 111)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 200, 5)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 210)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 201, 5)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 211)

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
        features, offsets = weightedHistogramBinFeatures(histogram[0], 'query',
                                                         TEST_WEIGHTS)
        self.assertEqual(set([queryLandmark1, queryTrigPoint1,
                              queryLandmark2, queryTrigPoint2]),
                         features)
        self.assertEqual({100: [1.5], 101: [1.5], 102: [1.5],
                          103: [1.5], 104: [1.5], 110: [1],
                          200: [1.5], 201: [1.5], 202: [1.5],
                          203: [1.5], 204: [1.5], 210: [1]}, offsets)

    def testWeightedOverlappingFeaturesDifferentWeights(self):
        """
        If two features overlap with different weights, the right offset dict
        must be returned.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 100, 5)
        queryTrigPoint1 = TrigPoint('Peaks', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 101, 5)
        subjectTrigPoint1 = TrigPoint('Peaks', 'P', 111)

        queryLandmark2 = Landmark('AlphaHelix_pi', 'A', 102, 5)
        queryTrigPoint2 = TrigPoint('Peaks', 'P', 112)
        subjectLandmark2 = Landmark('AlphaHelix_pi', 'A', 202, 5)
        subjectTrigPoint2 = TrigPoint('Peaks', 'P', 212)

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
        features, offsets = weightedHistogramBinFeatures(histogram[0], 'query',
                                                         TEST_WEIGHTS)
        self.assertEqual(set([queryLandmark1, queryTrigPoint1,
                              queryLandmark2, queryTrigPoint2]),
                         features)
        self.assertEqual({100: [1.5], 101: [1.5], 102: [1.5, 1],
                          103: [1.5, 1], 104: [1.5, 1], 105: [1], 106: [1],
                          110: [1], 112: [1]}, offsets)


class TestGetWeightedOffsets(TestCase):
    """
    Tests for the light.bin_score.getWeightedOffsets function.
    """
    def testEmptyDict(self):
        """
        An empty dict must return a count of 0.
        """
        result = getWeightedOffsets({})
        self.assertEqual(0, result)

    def testDictOneWeightPerKey(self):
        """
        An input dict with one weight per key must return the right count.
        """
        result = getWeightedOffsets({1: [1], 2: [1.5]})
        self.assertEqual(2.5, result)

    def testDictTwoWeightsPerKey(self):
        """
        An input dict which has two values for a key must return the right
        count.
        """
        result = getWeightedOffsets({1: [1], 2: [1, 1.5]})
        self.assertEqual(2.5, result)


class TestWeightedFeatureAAScore(TestCase):
    """
    Tests for the light.bin_score.WeightedFeatureAAScore class.
    """

    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of 0.0 if the query and
        subject both have no features.
        """
        histogram = Histogram()
        histogram.finalize()
        params = Parameters([], [])
        query = AARead('id1', 'A')
        subject = Subject('id2', 'A', 0)
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 0,
                'denominatorSubject': 0,
                'matchedOffsetCount': 0,
                'matchedQueryOffsetCount': 0,
                'weightedMatchedQueryOffsetCount': 0,
                'weightedMatchedSubjectOffsetCount': 0,
                'matchedRegionScore': 0.0,
                'matchedSubjectOffsetCount': 0,
                'maxQueryOffset': None,
                'maxSubjectOffset': None,
                'minQueryOffset': None,
                'minSubjectOffset': None,
                'numeratorQuery': 0,
                'numeratorSubject': 0,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 0,
            },
            analysis)

    def testEmptyBinQueryHasOneFeature(self):
        """
        A bin containing no hashes must have a score of zero, even if the query
        has a feature (but no hashes). There is no match region, so that part
        of the score is zero which causes the overall score to be zero.
        """
        histogram = Histogram(1)
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 0,
                'denominatorSubject': 0,
                'matchedOffsetCount': 0,
                'matchedQueryOffsetCount': 0,
                'weightedMatchedQueryOffsetCount': 0,
                'weightedMatchedSubjectOffsetCount': 0,
                'matchedRegionScore': 0.0,
                'matchedSubjectOffsetCount': 0,
                'maxQueryOffset': None,
                'maxSubjectOffset': None,
                'minQueryOffset': None,
                'minSubjectOffset': None,
                'numeratorQuery': 0,
                'numeratorSubject': 0,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 0,
            },
            analysis)

    def testEmptyBinQueryAndSubjectHaveOneFeature(self):
        """
        A bin containing no hashes must have a score of zero, when the query
        and subject both have one feature (but no hashes). There is no match
        region, so that part of the score is zero which causes the overall
        score to be zero.
        """
        histogram = Histogram(1)
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 'FRRRFRRRF')
        subject = Subject('id2', 'AAAAAAAAAAAAAAFRRRFRRRF', 0)
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(0.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 0,
                'denominatorSubject': 0,
                'matchedOffsetCount': 0,
                'matchedQueryOffsetCount': 0,
                'weightedMatchedQueryOffsetCount': 0,
                'weightedMatchedSubjectOffsetCount': 0,
                'matchedRegionScore': 0.0,
                'matchedSubjectOffsetCount': 0,
                'maxQueryOffset': None,
                'maxSubjectOffset': None,
                'minQueryOffset': None,
                'minSubjectOffset': None,
                'numeratorQuery': 0,
                'numeratorSubject': 0,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 0,
            },
            analysis)

    def testOneHashInBin(self):
        """
        A bin containing one hash must have a score of 1.0 if the query and
        subject have no additional (non-matching) hashes.
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      TEST_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 20,
                'matchedOffsetCount': 60.0,
                'matchedQueryOffsetCount': 20,
                'weightedMatchedQueryOffsetCount': 30.0,
                'weightedMatchedSubjectOffsetCount': 30.0,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 60.0,
            },
            analysis)

    def testOneHashInBinOccurringInTwoPlaces(self):
        """
        A bin containing one hash that occurs in two places must have a score
        of 1.0 if the query and subject have no additional (non-matching)
        hashes.
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      TEST_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 20,
                'matchedOffsetCount': 58.0,
                'matchedQueryOffsetCount': 20,
                'weightedMatchedSubjectOffsetCount': 29.0,
                'weightedMatchedQueryOffsetCount': 29.0,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 30,
                'maxSubjectOffset': 30,
                'minQueryOffset': 0,
                'minSubjectOffset': 1,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 58.0,
            },
            analysis)

    def testOneHashInBinQueryHasOneHashOutsideMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has an
        additional hash that is outside the match area (because the subject
        should be used to do the normalisation by length).
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      TEST_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 10,
                'matchedOffsetCount': 29.0,
                'matchedQueryOffsetCount': 10,
                'weightedMatchedSubjectOffsetCount': 14.5,
                'weightedMatchedQueryOffsetCount': 14.5,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 10,
                'maxQueryOffset': 10,
                'maxSubjectOffset': 10,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'numeratorQuery': 10,
                'numeratorSubject': 10,
                'normaliserQuery': 0.5,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 29.0,
            },
            analysis)

    def testOneHashInBinQueryHasTwoHashesOutsideMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has two
        additional hashes that are outside the match area (because the subject
        should be used to do the normalisation by length).
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 21,
                'denominatorSubject': 10,
                'matchedOffsetCount': 20.0,
                'matchedQueryOffsetCount': 10,
                'weightedMatchedSubjectOffsetCount': 10.0,
                'weightedMatchedQueryOffsetCount': 10.0,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 10,
                'maxQueryOffset': 110,
                'maxSubjectOffset': 110,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'numeratorQuery': 10,
                'numeratorSubject': 10,
                'normaliserQuery': 10 / 21,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 20.0,
            },
            analysis)

    def testOneHashInBinQuery2Subject1HashOutsideMatch(self):
        """
        A bin containing one hash must not have a score of 1.0 if the query has
        two and the subject one hash that are outside the match area.
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertAlmostEqual(2 / 3, score)
        self.assertEqual(
            {
                'denominatorQuery': 31,
                'denominatorSubject': 30,
                'matchedOffsetCount': 40.0,
                'matchedQueryOffsetCount': 20,
                'weightedMatchedSubjectOffsetCount': 20.0,
                'weightedMatchedQueryOffsetCount': 20.0,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 20 / 31,
                'normaliserSubject': 2 / 3,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 40.0,
            },
            analysis)

    def testOneHashInBinQuery1Subject2HashOutsideMatch(self):
        """
        A bin containing one hash must not have a score of 1.0 if the query has
        one and the subject two hashes that are outside the match area.
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertAlmostEqual(2 / 3, score)
        self.assertEqual(
            {
                'denominatorQuery': 30,
                'denominatorSubject': 31,
                'matchedOffsetCount': 40.0,
                'matchedQueryOffsetCount': 20,
                'weightedMatchedSubjectOffsetCount': 20.0,
                'weightedMatchedQueryOffsetCount': 20.0,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 119,
                'maxSubjectOffset': 119,
                'minQueryOffset': 100,
                'minSubjectOffset': 100,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 2 / 3,
                'normaliserSubject': 20 / 31,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 40.0,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedHashInsideMatch(self):
        """
        A bin containing one hash where the landmark and trig point do not
        overlap must have the correct score if the query has an additional
        non-matching hash that is inside the match area.
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual((21 + 21) / (21 + 21 + 10), score)
        self.assertEqual(
            {
                'denominatorQuery': 31,
                'denominatorSubject': 21,
                'matchedOffsetCount': 42.0,
                'matchedQueryOffsetCount': 21,
                'weightedMatchedSubjectOffsetCount': 21.0,
                'weightedMatchedQueryOffsetCount': 21.0,
                'matchedRegionScore': 42 / 52,
                'matchedSubjectOffsetCount': 21,
                'maxQueryOffset': 50,
                'maxSubjectOffset': 50,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'numeratorQuery': 31,
                'numeratorSubject': 21,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 52.0,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedHashExactlySpanningMatch(self):
        """
        A bin containing one hash must have a the correct score if the query
        has an additional hash that exactly spans the match area but the
        additional hashes' offsets match those of the match (and so do not
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual((9 + 9) / (9 + 9), score)
        self.assertEqual(
            {
                'denominatorQuery': 10,
                'denominatorSubject': 10,
                'matchedOffsetCount': 20.0,
                'matchedQueryOffsetCount': 10,
                'weightedMatchedSubjectOffsetCount': 10.0,
                'weightedMatchedQueryOffsetCount': 10.0,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 10,
                'maxQueryOffset': 13,
                'maxSubjectOffset': 13,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'numeratorQuery': 10,
                'numeratorSubject': 10,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 20.0,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedHashExceedingMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has
        an additional hash that exceeds the match area on both sides (because
        the subject is used for the score normalisation by length, the score
        is not affected).
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
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
                'weightedMatchedSubjectOffsetCount': 5,
                'weightedMatchedQueryOffsetCount': 5.0,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 5,
                'maxQueryOffset': 6,
                'maxSubjectOffset': 6,
                'minQueryOffset': 2,
                'minSubjectOffset': 2,
                'numeratorQuery': 5,
                'numeratorSubject': 5,
                'normaliserQuery': 0.5,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 10,
            },
            analysis)

    def testOneHashInBinQueryHasTwoUnmatchedFeaturesInsideMatch(self):
        """
        A bin containing one hash must have the correct score if the query has
        two additional features (making a hash) that are inside the match area
        but which do not overlap the features in the hash.
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(42 / 44, score)
        self.assertEqual(
            {
                'denominatorQuery': 23,
                'denominatorSubject': 21,
                'matchedOffsetCount': 42.0,
                'matchedQueryOffsetCount': 21,
                'weightedMatchedSubjectOffsetCount': 21.0,
                'weightedMatchedQueryOffsetCount': 21.0,
                'matchedRegionScore': 42 / 44,
                'matchedSubjectOffsetCount': 21,
                'maxQueryOffset': 50,
                'maxSubjectOffset': 50,
                'minQueryOffset': 0,
                'minSubjectOffset': 0,
                'numeratorQuery': 23,
                'numeratorSubject': 21,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 44.0,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchLeft(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has an
        additional feature that is only partly inside the match area (with the
        extra feature jutting out on the left).
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 25,
                'denominatorSubject': 20,
                'matchedOffsetCount': 40,
                'matchedQueryOffsetCount': 20,
                'weightedMatchedSubjectOffsetCount': 20,
                'weightedMatchedQueryOffsetCount': 20,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 20,
                'maxQueryOffset': 24,
                'maxSubjectOffset': 24,
                'minQueryOffset': 5,
                'minSubjectOffset': 5,
                'numeratorQuery': 20,
                'numeratorSubject': 20,
                'normaliserQuery': 0.8,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 40,
            },
            analysis)

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchRight(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has an
        additional feature that is only partly inside the match area (with the
        extra feature jutting out on the right).
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)
        self.assertEqual(
            {
                'denominatorQuery': 11,
                'denominatorSubject': 4,
                'matchedOffsetCount': 8,
                'matchedQueryOffsetCount': 4,
                'weightedMatchedSubjectOffsetCount': 4,
                'weightedMatchedQueryOffsetCount': 4,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 4,
                'maxQueryOffset': 5,
                'maxSubjectOffset': 5,
                'minQueryOffset': 2,
                'minSubjectOffset': 2,
                'numeratorQuery': 4,
                'numeratorSubject': 4,
                'normaliserQuery': 4 / 11,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 8,
            },
            analysis)

    def testTwoHashes(self):
        """
        A bin containing two hashes must have the correct score if the query
        and subject both have an additional feature inside their match areas.
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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        matched = (3 + 1) + (5 + 1) + (3 + 1) + (5 + 1)
        total = matched + (9 + 1) + (13 + 1)
        self.assertEqual(matched / total, score)
        self.assertEqual(
            {
                'denominatorQuery': 20,
                'denominatorSubject': 24,
                'matchedOffsetCount': 20.0,
                'matchedQueryOffsetCount': 10,
                'weightedMatchedSubjectOffsetCount': 10.0,
                'weightedMatchedQueryOffsetCount': 10.0,
                'matchedRegionScore': 20 / 44,
                'matchedSubjectOffsetCount': 10,
                'maxQueryOffset': 60,
                'maxSubjectOffset': 60,
                'minQueryOffset': 2,
                'minSubjectOffset': 2,
                'numeratorQuery': 20,
                'numeratorSubject': 24,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': score,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 44.0,
            },
            analysis)

    def testCompareEqualSequencesScoreMustBeOne(self):
        """
        If a sequence is compared to itself, the score must be 1.0. See
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
                                    scoreMethod='WeightedFeatureAAScore')
        result = db.find(pichninde, findParams, storeFullAnalysis=True)
        self.assertEqual(1.0, result.analysis[subjectIndex]['bestBinScore'])

        scoreAnalysis = result.analysis[
            subjectIndex]['significantBins'][0]['scoreAnalysis']
        self.maxDiff = None
        self.assertEqual(
            {
                'denominatorQuery': 210,
                'denominatorSubject': 210,
                'matchedOffsetCount': 420.0,
                'matchedQueryOffsetCount': 210,
                'weightedMatchedSubjectOffsetCount': 210.0,
                'weightedMatchedQueryOffsetCount': 210.0,
                'matchedRegionScore': 1.0,
                'matchedSubjectOffsetCount': 210,
                'maxQueryOffset': 290,
                'maxSubjectOffset': 290,
                'minQueryOffset': 1,
                'minSubjectOffset': 1,
                'numeratorQuery': 210,
                'numeratorSubject': 210,
                'normaliserQuery': 1.0,
                'normaliserSubject': 1.0,
                'score': 1.0,
                'scoreClass': WeightedFeatureAAScore,
                'totalOffsetCount': 420.0,
            },
            scoreAnalysis)

    def testScoresMustBeSymmetric(self):
        """
        When comparing two sequences, the scores must be the same, no matter
        which one is used as the query or subject.

        This was a (formerly) failing test built during the resolution of
        https://github.com/acorg/light-matter/issues/341 based on two of the
        sequences received from Sandra Junglen on March 13, 2015.
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
                                    scoreMethod='WeightedFeatureAAScore')

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
        faas = WeightedFeatureAAScore(histogram, query, subject, params,
                                      DEFAULT_WEIGHTS)
        score, analysis = faas.calculateScore(0)
        self.assertEqual(1.0, score)

        expected = (
            'Score method: WeightedFeatureAAScore\n'
            'Matched offset range in query: 5 to 24\n'
            'Matched offset range in subject: 5 to 24\n'
            'Total (query+subject) AA offsets in matched hashes: 40\n'
            'Subject AA offsets in matched hashes: 20\n'
            'Query AA offsets in matched hashes: 20\n'
            'Total (query+subject) AA offsets in hashes in matched '
            'region: 40\n'
            'Weighted Subject AA offsets in matched hashes: 20\n'
            'Weighted Query AA offsets in matched hashes: 20\n'
            'Matched region score 1.0000 (40 / 40)\n'
            'Query normalizer: 0.8000 (20 / 25)\n'
            'Subject normalizer: 1.0000 (20 / 20)\n'
            'Score: 1.0000')
        self.assertEqual(expected,
                         WeightedFeatureAAScore.printAnalysis(analysis))
        self.assertEqual(expected,
                         analysis['scoreClass'].printAnalysis(analysis))
