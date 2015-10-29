import warnings
from unittest import TestCase

from dark.reads import AARead

from light.features import Landmark, TrigPoint
from light.parameters import Parameters, FindParameters
from light.subject import Subject
from light.score import (
    MinHashesScore, FeatureMatchingScore, FeatureAAScore, histogramBinFeatures,
    featureInRange)
from light.histogram import Histogram
from light.landmarks import AlphaHelix


class TestMinHashesScore(TestCase):
    """
    Tests for the light.score.MinHashesScore class.
    """
    def testEmptyBin(self):
        """
        A bin containing no elements must have a score of 0.0
        """
        histogram = Histogram()
        histogram.finalize()
        mhs = MinHashesScore(histogram, 0)
        self.assertEqual(0.0, mhs.calculateScore(0))

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
        self.assertEqual(0.0, mhs.calculateScore(0))

    def testBinWithOneElementMinHashCountOne(self):
        """
        If the minimum hash count for the query/subject is one and the bin has
        one thing in it, the score for the bin must be 1.0
        """
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 1)
        self.assertEqual(1.0, mhs.calculateScore(0))

    def testBinWithOneElementMinHashCountTwo(self):
        """
        If the minimum hash count for the query/subject is two and the bin has
        one thing in it, the score for the bin must be 0.5
        """
        histogram = Histogram(1)
        histogram.add(555)
        histogram.finalize()
        mhs = MinHashesScore(histogram, 2)
        self.assertEqual(0.5, mhs.calculateScore(0))

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
            mhs.calculateScore(0)
            self.assertEqual(1, len(w))
            self.assertTrue(issubclass(w[0].category, RuntimeWarning))
            error = ('Bin contains 2 deltas for a query/subject pair with a '
                     'minimum hash count of only 1.')
            self.assertIn(error, str(w[0].message))


class TestHistogramBinFeatures(TestCase):
    """
    Tests for the light.score.histogramBinFeatures function.
    """
    def testInvalidQueryOrSubjectSpecifier(self):
        """
        If an invalid value for queryOrSubject ('xxx') is passed to
        histogramBinFeatures, it must raise a KeyError (when looking
        in the bin item dictionary for 'xxxLandmark').
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 111)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        self.assertRaisesRegexp(KeyError, 'xxxLandmark', histogramBinFeatures,
                                histogram[0], 'xxx')

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
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 111)
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
        self.assertEqual(set([100, 110, 119]), offsets)

    def testOneFeatureSubject(self):
        """
        If a histogram bin has just one feature, histogramBinFeatures must
        return the details of that feature in the subject.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 111)
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
        self.assertEqual(set([101, 111, 120]), offsets)

    def testOneFeatureTwoLocations(self):
        """
        If a histogram bin has one feature that appears in two places,
        histogramBinFeatures must return the details of that feature.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint1 = TrigPoint('Peak', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint1 = TrigPoint('Peak', 'P', 111)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        queryTrigPoint2 = TrigPoint('Peak', 'P', 210)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 201, 20)
        subjectTrigPoint2 = TrigPoint('Peak', 'P', 211)

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
        self.assertEqual(set([100, 110, 119, 200, 210, 219]), offsets)


class TestFeatureInRange(TestCase):
    """
    Tests for the light.score.featureInRange function.
    """
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


class TestFeatureMatchingScore(TestCase):
    """
    Tests for the light.score.FeatureMatchingScore class.
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
        self.assertEqual(0.0, fms.calculateScore(0))

    def testOneHashInBin(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four) if the query and subject have no
        additional (non-matching) features. The reward is multiplied by four
        because both the landmark and the trig point (2) appear in both the
        query and the subject (2) and 2 x 2 = 4.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 110)
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
        self.assertEqual(4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

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
        queryTrigPoint1 = TrigPoint('Peak', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint1 = TrigPoint('Peak', 'P', 111)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        queryTrigPoint2 = TrigPoint('Peak', 'P', 210)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 201, 20)
        subjectTrigPoint2 = TrigPoint('Peak', 'P', 211)

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
        self.assertEqual(8 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

    def testOneHashInBinQueryHasOneFeatureOutsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four if the query has an additional feature
        that is outside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 110)
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
        self.assertEqual(4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

    def testOneHashInBinQueryHasTwoFeaturesOutsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four if the query has two additional
        features that are outside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 110)
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
        self.assertEqual(4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

    def testOneHashInBinQueryAndSubjectHaveOneFeaturesOutsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature
        match reward multiplied by four if the query and the subject each have
        one additional feature that are both outside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 110)
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
        self.assertEqual(4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
                         fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureInsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus the feature mismatch score if the
        query has an additional feature that is inside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 10)
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
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE +
            FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureExactlySpanningMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus the feature mismatch score if the
        query has an additional feature that exactly spans the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        queryTrigPoint = TrigPoint('Peak', 'P', 5)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        subjectTrigPoint = TrigPoint('Peak', 'P', 5)
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
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE +
            FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureExceedingMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four if the query has an additional feature that
        exceeds the match area on both sides.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 2, 5)
        queryTrigPoint = TrigPoint('Peak', 'P', 5)
        subjectLandmark = Landmark('AlphaHelix', 'A', 2, 5)
        subjectTrigPoint = TrigPoint('Peak', 'P', 5)
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
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasTwoUnmatchedFeatureInsideMatch(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus twice the feature mismatch score if
        the query has two additional features that are inside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 10)
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
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE +
            2 * FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchLeft(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four if the query has an additional feature that
        is only partly inside the match area (with the extra feature jutting
        out on the left).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 10)
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
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchRight(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four if the query has an additional feature that
        is only partly inside the match area (with the extra feature jutting
        out on the right).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 2, 3)
        queryTrigPoint = TrigPoint('Peak', 'P', 5)
        subjectLandmark = Landmark('AlphaHelix', 'A', 2, 3)
        subjectTrigPoint = TrigPoint('Peak', 'P', 5)
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
        self.assertEqual(
            4 * FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
            fms.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureInMatchNonDefault(self):
        """
        A bin containing one hash must have a score that is the feature match
        reward multiplied by four, minus the feature mismatch score if the
        query has an additional feature that is inside the match area,
        including when non-default values are used for the feature match and
        mismatch scores.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 10)
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
        self.assertEqual(4 * 3.1 - 1.2, fms.calculateScore(0))


class TestFeatureAAScore(TestCase):
    """
    Tests for the light.score.FeatureAAScore class.
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
        self.assertEqual(0.0, faas.calculateScore(0))

    def testOneHashInBin(self):
        """
        A bin containing one hash must have a score of 1.0 if the query and
        subject have no additional (non-matching) features.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 110)
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
        self.assertEqual(1.0, faas.calculateScore(0))

    def testOneHashInBinOccurringInTwoPlaces(self):
        """
        A bin containing one hash that occurs in two place must have a score
        of 1.0 if the query and subject have no additional (non-matching)
        features.
        """
        queryLandmark1 = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint1 = TrigPoint('Peak', 'P', 110)
        subjectLandmark1 = Landmark('AlphaHelix', 'A', 101, 20)
        subjectTrigPoint1 = TrigPoint('Peak', 'P', 111)

        queryLandmark2 = Landmark('AlphaHelix', 'A', 200, 20)
        queryTrigPoint2 = TrigPoint('Peak', 'P', 210)
        subjectLandmark2 = Landmark('AlphaHelix', 'A', 201, 20)
        subjectTrigPoint2 = TrigPoint('Peak', 'P', 211)

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
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual(1.0, faas.calculateScore(0))

    def testOneHashInBinQueryHasOneFeatureOutsideMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has an
        additional feature that is outside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 110)
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
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual(1.0, faas.calculateScore(0))

    def testOneHashInBinQueryHasTwoFeaturesOutsideMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has two
        additional features that are outside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 110)
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
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual(1.0, faas.calculateScore(0))

    def testOneHashInBinQueryAndSubjectHaveOneFeaturesOutsideMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query and the
        subject each have one additional feature that are both outside the
        match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 110)
        subjectLandmark = Landmark('AlphaHelix', 'A', 100, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 110)
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
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual(1.0, faas.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureInsideMatch(self):
        """
        A bin containing one hash where the landmark and trig point do not
        overlap must have the correct score if the query has an additional
        feature that is inside the match area.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 50)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 50)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 20 * 'A' + 'FRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual((21 + 21) / (21 + 21 + 9),
                         faas.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureExactlySpanningMatch(self):
        """
        A bin containing one hash must have a the correct score if the query
        has an additional feature that exactly spans the match area but the
        additional feature's offsets match those of the match (and so do not
        affect the score).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        queryTrigPoint = TrigPoint('Peak', 'P', 5)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 9)
        subjectTrigPoint = TrigPoint('Peak', 'P', 5)
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
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual(
            (9 + 9) / (9 + 9),
            faas.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureExceedingMatch(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has
        an additional feature that exceeds the match area on both sides.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 2, 5)
        queryTrigPoint = TrigPoint('Peak', 'P', 5)
        subjectLandmark = Landmark('AlphaHelix', 'A', 2, 5)
        subjectTrigPoint = TrigPoint('Peak', 'P', 5)
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
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual(1.0, faas.calculateScore(0))

    def testOneHashInBinQueryHasTwoUnmatchedFeatureInsideMatch(self):
        """
        A bin containing one hash must have the correct score if the query has
        two additional features that are inside the match area but which do not
        overlap the features in the hash.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 50)
        subjectLandmark = Landmark('AlphaHelix', 'A', 0, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 50)
        histogram = Histogram(1)
        histogram.add(44, {
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        })
        histogram.finalize()
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 22 * 'A' + 'FRRRFRRRFAFRRRFRRRF')
        subject = Subject('id2', 'A', 0)
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual(
            ((20 + 1) + (20 + 1)) / ((20 + 1) + (20 + 1) + (2 * 9)),
            faas.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchLeft(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has an
        additional feature that is only partly inside the match area (with the
        extra feature jutting out on the left).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        queryTrigPoint = TrigPoint('Peak', 'P', 10)
        subjectLandmark = Landmark('AlphaHelix', 'A', 5, 20)
        subjectTrigPoint = TrigPoint('Peak', 'P', 10)
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
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual(1.0, faas.calculateScore(0))

    def testOneHashInBinQueryHasOneUnmatchedFeatureOverlappingMatchRight(self):
        """
        A bin containing one hash must have a score of 1.0 if the query has an
        additional feature that is only partly inside the match area (with the
        extra feature jutting out on the right).
        """
        queryLandmark = Landmark('AlphaHelix', 'A', 2, 3)
        queryTrigPoint = TrigPoint('Peak', 'P', 5)
        subjectLandmark = Landmark('AlphaHelix', 'A', 2, 3)
        subjectTrigPoint = TrigPoint('Peak', 'P', 5)
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
        faas = FeatureAAScore(histogram, query, subject, params)
        self.assertEqual(1.0, faas.calculateScore(0))

    def testTwoHashes(self):
        """
        A bin containing two hashes must have the correct score if the query
        and subject both have an additional feature inside their match areas.
        """
        queryLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        queryTrigPoint1 = TrigPoint('Peak', 'P', 10)
        subjectLandmark1 = Landmark('AlphaHelix_pi', 'C', 2, 3)
        subjectTrigPoint1 = TrigPoint('Peak', 'P', 10)

        queryLandmark2 = Landmark('AlphaHelix_pi', 'C', 50, 5)
        queryTrigPoint2 = TrigPoint('Peak', 'P', 60)
        subjectLandmark2 = Landmark('AlphaHelix_pi', 'C', 50, 5)
        subjectTrigPoint2 = TrigPoint('Peak', 'P', 60)

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
        params = Parameters([AlphaHelix], [])
        query = AARead('id', 20 * 'A' + 'FRRRFRRRF')
        subject = Subject('id2', 25 * 'A' + 'FRRRFRRRFRRRF', 0)
        faas = FeatureAAScore(histogram, query, subject, params)
        matched = (3 + 1) + (5 + 1) + (3 + 1) + (5 + 1)
        total = matched + 9 + 13
        self.assertEqual(matched / total, faas.calculateScore(0))
